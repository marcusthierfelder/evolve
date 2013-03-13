package evolve

import (
	"fmt"
	"log"
	"math"
	"reflect"
	"strconv"
)

var (
	tmp_vl   []VarList
	time_int = "none"

	rk_eps_min = 1e-4
	rk_eps_max = 1e-3
	rk_it_max  = 20

	rhs_ptr interface{}
)

type Grid struct {
	dim              int     // dimensions used
	nx, ny, nz       int     // number of points
	ntot, di, dj, dk int     // number of points
	gh               int     // number of ghosts
	time             float64 // time ...
	dx, dy, dz       float64 // grid spacing
	x0, y0, z0       float64 // coordinate offset
	field            []field // data storage

	boundary, boundaryA, boundaryB, boundaryC []int
	outvars                        []string
}
type field struct {
	name  string
	sync  bool
	data  []float64
	gtype byte
}
type VarList struct {
	field []*field
	grid  *Grid
}

// create a grid
func CreateGrid(Nx []int, x0, dx []float64) *Grid {
	var grid Grid

	grid.dim = len(Nx)
	if grid.dim < 0 || grid.dim > 3 {
		log.Fatal("use only 0D 1D 2D 3D as grid size")
	}
	if len(Nx) != len(x0) || len(Nx) != len(dx) {
		log.Fatal("wrong dimensions used")
	}

	tNx := []int{1, 1, 1}
	tdx := []float64{1., 1., 1.}
	tx0 := []float64{0., 0., 0.}

	copy(tNx, Nx)
	copy(tdx, dx)
	copy(tx0, x0)

	grid.nx, grid.ny, grid.nz = tNx[0], tNx[1], tNx[2]
	grid.dx, grid.dy, grid.dz = tdx[0], tdx[1], tdx[2]
	grid.x0, grid.y0, grid.z0 = tx0[0], tx0[1], tx0[2]

	grid.ntot = tNx[0] * tNx[1] * tNx[2]
	grid.di = 1
	grid.dj = tNx[0]
	grid.dk = tNx[0] * tNx[1]
	grid.time = 0.

	grid.AddVars([]string{"x", "y", "z"})
	x := grid.GetVar("x")
	y := grid.GetVar("y")
	z := grid.GetVar("z")
	ijk := 0
	for k := 0; k < tNx[2]; k++ {
		for j := 0; j < tNx[1]; j++ {
			for i := 0; i < tNx[0]; i++ {
				x[ijk] = tx0[0] + tdx[0]*float64(i)
				y[ijk] = tx0[1] + tdx[1]*float64(j)
				z[ijk] = tx0[2] + tdx[2]*float64(k)

				ijk++
			}
		}
	}

	grid.boundary = make([]int, grid.ntot)
	grid.boundaryA = make([]int, (grid.nx+1)*grid.ny*grid.nz)
	grid.boundaryB = make([]int, grid.nx*(grid.ny+1)*grid.nz)
	grid.boundaryC = make([]int, grid.nx*grid.ny*(grid.nz+1))

	return &grid
}

// get back the size of the grid 
// Nx number of pts in one direction
// di offset to get to the next point in one direction within the global stack
// dx gridspacing in one direction
func (grid *Grid) GetSize() (Nx,di []int, dx []float64) {

	Nx = []int{grid.nx, grid.ny, grid.nz}
	di = []int{grid.di, grid.dj, grid.dj}
	dx = []float64{grid.dx, grid.dy, grid.dz}

	return Nx[0:grid.dim], di[0:grid.dim], dx[0:grid.dim]
}

// get back the current time
func (grid *Grid) GetTime() float64 {
	return grid.time
}

// time integrator initialization which have to be called before
// the actual integration in order to set integrator and pass rhs
// computation function 
func (grid *Grid) TimeInt_init(uc VarList, integrator string, rhs_func interface{}) {
	// better too much ... its the simplest way 
	tmp_vl = []VarList{
		grid.CreateVarlist(), grid.CreateVarlist(), grid.CreateVarlist(),
		grid.CreateVarlist(), grid.CreateVarlist(), grid.CreateVarlist(),
		grid.CreateVarlist(), grid.CreateVarlist(), grid.CreateVarlist()}

	for _, v := range uc.field {
		switch integrator {
		case "euler":
			grid.AddVar(v.name + "_r")
			tmp_vl[0].AddVar(v.name + "_r")
		case "rk4", "icn":
			grid.AddVar(v.name + "_p")
			tmp_vl[0].AddVar(v.name + "_p")
			grid.AddVar(v.name + "_r")
			tmp_vl[1].AddVar(v.name + "_r")
			grid.AddVar(v.name + "_v")
			tmp_vl[2].AddVar(v.name + "_v")
		case "rk45":
			grid.AddVar(v.name + "_v")
			tmp_vl[0].AddVar(v.name + "_v")
			for k := 1; k <= 6; k++ {
				grid.AddVar(v.name + "_k" + strconv.Itoa(k))
				tmp_vl[k].AddVar(v.name + "_k" + strconv.Itoa(k))
			}
			grid.AddVar(v.name + "_rk4")
			tmp_vl[7].AddVar(v.name + "_rk4")
			grid.AddVar(v.name + "_rk5")
			tmp_vl[8].AddVar(v.name + "_rk5")
		default:
			log.Fatal("use euler icl rk4 rk45 ... or implement more ;)")

		}
	}

	// set integrator
	time_int = integrator
	rhs_ptr = rhs_func
}

// set the rk4 accuracy options 
func (grid *Grid) SetRK45Accuracy(min, max float64, it int) {
	rk_eps_min = min
	rk_eps_max = max
	rk_it_max = it
}

// actual timeintegrator, which does (hopefully) everything by itself
func (grid *Grid) TimeInt(uc VarList, dt *float64) {

	switch time_int {
	case "euler":
		grid.euler(uc, *dt)
	case "icn":
		grid.icn(uc, *dt)
	case "rk4":
		grid.rk4(uc, *dt)
	case "rk45":
		*dt = grid.rk45(uc, *dt)
	default:
		panic("integrator not implemented")
	}

	grid.time += *dt
}

/* helper functions for the time integrator */
func (grid *Grid) cpy(vl1 VarList, vl2 VarList) {
	for i, v1 := range vl1.field {
		v2 := vl2.field[i]

		ntot := len(v1.data)
		for ijk := 0; ijk < ntot; ijk++ {
			v1.data[ijk] = v2.data[ijk]
		}
	}
}

func (grid *Grid) addto(vl1 VarList, c float64, vl2 VarList) {
	for i, v1 := range vl1.field {
		v2 := vl2.field[i]

		ntot := len(v1.data)
		for ijk := 0; ijk < ntot; ijk++ {
			v1.data[ijk] += c * v2.data[ijk]
		}
	}
}

func (grid *Grid) add(vl1 VarList, c2 float64, vl2 VarList, c3 float64, vl3 VarList) {
	for i, v1 := range vl1.field {
		v2 := vl2.field[i]
		v3 := vl3.field[i]

		ntot := len(v1.data)
		for ijk := 0; ijk < ntot; ijk++ {
			v1.data[ijk] = c2*v2.data[ijk] + c3*v3.data[ijk]
		}
	}
}

func (grid *Grid) vlnorm(vl1 VarList, vl2 VarList) float64 {

	m := 0.
	for i, v1 := range vl1.field {
		v2 := vl2.field[i]

		min := math.Abs(v1.data[0] - v2.data[0])
		max := 0.

		ntot := len(v1.data)
		for ijk := 0; ijk < ntot; ijk++ {
			v := math.Abs(v1.data[ijk] - v2.data[ijk])
			if v < min {
				min = v
			}
			if v > max {
				max = v
			}
		}

		//fmt.Println(min,max)
		if max > m {
			m = max
		}
	}

	return m
}

func (grid *Grid) rhs(r, uc VarList) {
	f := reflect.ValueOf(rhs_ptr)
	in := make([]reflect.Value, 3)
	in[0] = reflect.ValueOf(grid)
	in[1] = reflect.ValueOf(r)
	in[2] = reflect.ValueOf(uc)
	f.Call(in)
}

/* time integratos */
func (grid *Grid) rk4(uc VarList, dt float64) {

	p := tmp_vl[0]
	r := tmp_vl[1]
	v := tmp_vl[2]

	grid.cpy(p, uc)

	grid.rhs(r, uc)
	grid.addto(uc, dt/6., r)

	grid.add(v, 1.0, p, dt/2., r)
	grid.rhs(r, v)
	grid.addto(uc, dt/3., r)

	grid.add(v, 1.0, p, dt/2., r)
	grid.rhs(r, v)
	grid.addto(uc, dt/3., r)

	grid.add(v, 1.0, p, dt, r)
	grid.rhs(r, v)
	grid.addto(uc, dt/6., r)
}

func (grid *Grid) rk45(uc VarList, dt float64) float64 {

	v := tmp_vl[0]
	k1 := tmp_vl[1]
	k2 := tmp_vl[2]
	k3 := tmp_vl[3]
	k4 := tmp_vl[4]
	k5 := tmp_vl[5]
	k6 := tmp_vl[6]
	rk4 := tmp_vl[7]
	rk5 := tmp_vl[8]

	for k := 0; k < rk_it_max; k++ {

		grid.rhs(k1, uc)

		grid.cpy(v, uc)
		grid.addto(v, 1./4.*dt, k1)
		grid.rhs(k2, v)

		grid.cpy(v, uc)
		grid.addto(v, 2./32.*dt, k1)
		grid.addto(v, 9./32.*dt, k2)
		grid.rhs(k3, v)

		grid.cpy(v, uc)
		grid.addto(v, 1932./2197.*dt, k1)
		grid.addto(v, -7200./2197.*dt, k2)
		grid.addto(v, 7296./2197.*dt, k3)
		grid.rhs(k4, v)

		grid.cpy(v, uc)
		grid.addto(v, 439./216.*dt, k1)
		grid.addto(v, -8./1.*dt, k2)
		grid.addto(v, 3680./513.*dt, k3)
		grid.addto(v, -845./4104.*dt, k4)
		grid.rhs(k5, v)

		grid.cpy(v, uc)
		grid.addto(v, -8./27.*dt, k1)
		grid.addto(v, 2./1.*dt, k2)
		grid.addto(v, -3544./2565.*dt, k3)
		grid.addto(v, 1859./4104.*dt, k4)
		grid.addto(v, -11./40.*dt, k5)
		grid.rhs(k6, v)

		// RK4
		grid.cpy(rk4, uc)
		grid.addto(rk4, 25./216.*dt, k1)
		grid.addto(rk4, 1408./2565.*dt, k3)
		grid.addto(rk4, 2197./4104.*dt, k4)
		grid.addto(rk4, -1./5.*dt, k5)

		// RK5
		grid.cpy(rk5, uc)
		grid.addto(rk5, 16./135.*dt, k1)
		grid.addto(rk5, 6656./12825.*dt, k3)
		grid.addto(rk5, 28561./56430.*dt, k4)
		grid.addto(rk5, -18./100.*dt, k5)
		grid.addto(rk5, -2./55.*dt, k6)

		eps := grid.vlnorm(rk4, rk5)
		//break

		if eps < rk_eps_min {
			dt *= 2.
		} else if eps > rk_eps_max {
			dt *= 0.5
		} else {
			break
		}
		fmt.Println("changed timestep to ", dt, "epsl= ", eps)

	}

	grid.cpy(uc, rk5)

	return dt
}

func (grid *Grid) icn(uc VarList, dt float64) {

	p := tmp_vl[0]
	r := tmp_vl[1]
	v := tmp_vl[2]

	grid.cpy(p, uc)

	for i := 0; i < 3; i++ {
		grid.add(v, 0.5, uc, 0.5, p)
		grid.rhs(r, v)
		grid.add(uc, 1.0, p, dt, r)
	}
}

func (grid *Grid) euler(uc VarList, dt float64) {

	r := tmp_vl[0]

	grid.rhs(r, uc)
	grid.addto(uc, dt, r)
}

// add a variable to global stuck
func (grid *Grid) AddVar(name string) {
	fmt.Println("AddVar: ", name)

	if false {
		fmt.Println(grid.field, len(grid.field))
	}

	l := len(grid.field)
	tmp := make([]field, l+1)
	copy(tmp, grid.field)
	grid.field = tmp

	f := field{
		name: name,
		data: make([]float64, grid.nx*grid.ny)}

	grid.field[l] = f
}

// add a variables to global stuck
func (grid *Grid) AddVars(names []string) {
	for _, name := range names {
		grid.AddVar(name)
	}
}

// get the data of this field
func (grid *Grid) GetVar(name string) []float64 {
	ptr := grid.getfield(name)
	return ptr.data
}

// get the boundary data as a screen
func (grid *Grid) GetBound(name string) []int {
	switch name {
	case "":
		return grid.boundary
	case "A":
		return grid.boundaryA
	case "B":
		return grid.boundaryB
	case "C":
		return grid.boundaryC
	default:
			panic("boundary does not exist!")
	}

	return []int{}
}

func (grid *Grid) getfield(name string) *field {
	var ptr *field

	ptr = nil
	for i, f := range grid.field {
		if f.name == name {
			ptr = &(grid.field[i])
			// this DOES NOT WORK, why?
			// ptr =&f
		}
	}

	if ptr == nil {
		log.Fatal("var \"" + name + "\" does not exist")
	}
	return ptr
}

//manage variables which have to be passed to the timeintegrators 
func (grid *Grid) CreateVarlist() VarList {
	var vl VarList
	vl.grid = grid

	return vl
}

// add a variable to this varlist
func (vl *VarList) AddVar(name string) {
	l := len(vl.field)
	tmp := make([]*field, l+1)
	copy(tmp, vl.field)
	vl.field = tmp

	vl.field[l] = vl.grid.getfield(name)
}

// add a variables to this varlist
func (vl *VarList) AddVars(names []string) {
	for _, name := range names {
		vl.AddVar(name)
	}
}

// get the data of the i-th stored variable in this list
func (vl *VarList) GetVar(i int) []float64 {
	if i >= len(vl.field) {
		log.Fatal("GetVar, number out of list")
	}
	return vl.field[i].data
}
