package evolve

import (
  _ "flag"
	"fmt"
	"io"
	"log"
	"math"
	"os"
	 "reflect"
	"strconv"
)




var (
	tmp_vl     []VarList
	INTEGRATOR = "none"

	rk_eps_min = 1e-4
	rk_eps_max = 1e-3
	rk_it_max  = 20

	rhs_ptr interface{}

)


type Grid struct {
	nx, ny, nz   int // ... as in 
	ntot, di, dj, dk int
	gh           int // number of ghosts
	nt, ot       int
	time         float64 // time ...
	dx, dy, dz       float64
x0, y0, z0       float64

	field                          []Field // data storage
	boundary, boundaryA, boundaryB []int

	outvars []string
}
type Field struct {
	name  string
	sync  bool
	data  []float64
	gtype byte
}
type VarList struct {
	field []*Field
	grid  *Grid
}


/* set Gridvalue */
func CreateGrid(Nx []int, x0,dx []float64) (*Grid) {
  	var grid Grid

  	grid.nx = Nx[0]
grid.x0 = x0[0]
grid.dx = dx[0]
grid.ntot = Nx[0]

grid.nt = 10
grid.ot = 1


  	return &grid
}


/* time integrator wrapers */

func (grid *Grid)  Testcall(m map[string]interface{}, name string) {
	f := reflect.ValueOf(m[name])


	fmt.Println(f)
	fmt.Println(reflect.TypeOf(f))
	fmt.Println(reflect.TypeOf(m[name]))

	in := make([]reflect.Value, 0)
	f.Call(in)

}


func (grid *Grid)  Testcall2(i interface{}) {
	f := reflect.ValueOf(i)

	in := make([]reflect.Value, 0)
	f.Call(in)

	rhs_ptr = i
}

func (grid *Grid)  Testcall3() {
	f := reflect.ValueOf(rhs_ptr)
	in := make([]reflect.Value, 0)
	fmt.Println("hello---")
	f.Call(in)
}






func (grid *Grid) rhs(r, evl VarList) {
	f := reflect.ValueOf(rhs_ptr)
	in := make([]reflect.Value, 2)
	in[0] = reflect.ValueOf(r)
	in[1] = reflect.ValueOf(evl)
	f.Call(in)
}

func (grid *Grid) TimeInt_init(uc VarList, integrator string) {
	tmp_vl = []VarList{
		grid.vlalloc(), grid.vlalloc(), grid.vlalloc(),
		grid.vlalloc(), grid.vlalloc(), grid.vlalloc(),
		grid.vlalloc(), grid.vlalloc(), grid.vlalloc()}

	for _, v := range uc.field {
		switch integrator {
		case "euler":
			grid.AddVar(v.name + "_r")
			tmp_vl[1].AddVar(v.name + "_r")
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
		}
	}

	INTEGRATOR = integrator
}

func (grid *Grid) TimeInt(uc VarList, dt *float64) {

	switch INTEGRATOR {
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

/* variable storage stuff */
func (grid *Grid) AddVar(name string) {
	fmt.Println("AddVar: ", name)

	if false {
		fmt.Println(grid.field, len(grid.field))
	}

	l := len(grid.field)
	tmp := make([]Field, l+1)
	copy(tmp, grid.field)
	grid.field = tmp

	f := Field{
		name: name,
		data: make([]float64, grid.nx*grid.ny)}

	grid.field[l] = f
}

func (grid *Grid) AddVars(names []string) {
	for _, name := range names {
		grid.AddVar(name)
	}
}

func (grid *Grid) GetVar(name string) []float64 {
	ptr := grid.GetField(name)
	return ptr.data
}

func (grid *Grid) GetField(name string) *Field {
	var ptr *Field

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

/* varlist stuff */
func (grid *Grid) vlalloc() VarList {
	var vl VarList
	vl.grid = grid

	return vl
}

func (vl *VarList) AddVar(name string) {
	l := len(vl.field)
	tmp := make([]*Field, l+1)
	copy(tmp, vl.field)
	vl.field = tmp

	vl.field[l] = vl.grid.GetField(name)
}

func (vl *VarList) AddVars(names []string) {
	for _, name := range names {
		vl.AddVar(name)
	}
}

func (vl *VarList) GetVar(i int) []float64 {
	if i >= len(vl.field) {
		log.Fatal("GetVar, number out of list")
	}
	return vl.field[i].data
}

/* output */
func (grid *Grid) output_vtk(data string, file string, it int) {

	/* init buffer */
	//x := grid.GetVar("x")
	//y := grid.GetVar("y")
	d := grid.GetVar(data)

	flag := os.O_CREATE | os.O_TRUNC | os.O_RDWR

	file = fmt.Sprintf("%s_%04d.vtk", file, it)

	f, err := os.OpenFile(file, flag, 0666)
	if err != nil {
		log.Fatal(err)
		return
	}

	fmt.Println("write ", file, grid.time)
	//io.WriteString(f, fmt.Sprintf("#Time = %e\n", grid.time))
	io.WriteString(f, fmt.Sprintf("# vtk DataFile Version 2.0\n"))
	io.WriteString(f, fmt.Sprintf("variable f, level 0, time %16.9e\n", grid.time))
	io.WriteString(f, fmt.Sprintf("ASCII\n"))
	io.WriteString(f, fmt.Sprintf("DATASET STRUCTURED_POINTS\n"))
	io.WriteString(f, fmt.Sprintf("DIMENSIONS %d %d %d\n", grid.nx, grid.ny, grid.nz))
	io.WriteString(f, fmt.Sprintf("ORIGIN  %16.9e %16.9e %16.9e\n", 0., 0., 0.))
	io.WriteString(f, fmt.Sprintf("SPACING %16.9e %16.9e %16.9e\n", grid.dx, grid.dy, grid.dz))
	io.WriteString(f, fmt.Sprintf("\n"))
	io.WriteString(f, fmt.Sprintf("POINT_DATA %d\n", grid.nx*grid.ny))
	io.WriteString(f, fmt.Sprintf("SCALARS scalars float\n"))
	io.WriteString(f, fmt.Sprintf("LOOKUP_TABLE default\n"))

	ij := 0
	for k := 0; k < grid.nz; k++ {
	for j := 0; j < grid.ny; j++ {
		for i := 0; i < grid.nx; i++ {
			io.WriteString(f, fmt.Sprintf("%16.09e\n", d[ij]))
			ij++
		}
	}
}
	//io.WriteString(f, "\n")

	f.Close()
}