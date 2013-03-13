package evolve

import (
	"code.google.com/p/lvd.go/cdf"
	"fmt"
	"os"
)

var (
	outfile_netcdf string
)

func (grid *Grid) Output_netcdf_init(nout int) {

	if grid.dim < 1 || grid.dim > 3 {
		panic("this dimension is not supported")
	}

	outfile_netcdf = outpath + "/" + outfile + ".ncf"
	fmt.Println("Init netcdf file:", outfile_netcdf)

	// could be done better, however this f***ing permutation of the dimensions 
	// have to be considered
	var h *cdf.Header
	var c1, c2 []string
	var d1, d2 []int
	switch grid.dim {
	case 1:
		c1 = []string{"x"}
		c2 = []string{"time", "x"}
		d1 = []int{nout, grid.nx}
		d2 = []int{0}
	case 2:
		c1 = []string{"y", "x"}
		c2 = []string{"time", "y", "x"}
		d1 = []int{nout, grid.ny, grid.nx}
		d2 = []int{0, 0}
	case 3:
		c1 = []string{"z", "y", "x"}
		c2 = []string{"time", "z", "y", "x"}
		d1 = []int{nout, grid.nz, grid.ny, grid.nx}
		d2 = []int{0, 0, 0}
	}

	// generic ndim header
	h = cdf.NewHeader(c2, d1)
	h.AddVariable("x", []string{"x"}, []float32{0.})
	if grid.dim > 1 {
		h.AddVariable("y", []string{"y"}, []float32{0.})
	}
	if grid.dim > 2 {
		h.AddVariable("z", []string{"z"}, []float32{0.})
	}
	h.AddVariable("boundary", c1, []int32{0})
	for _, vname := range grid.outvars {
		h.AddVariable(vname, c2, []float32{0.})
	}

	h.Define()
	ff, _ := os.Create(outfile_netcdf)
	f, _ := cdf.Create(ff, h)

	// things which are independent 
	buf := make([]int32, grid.ntot)
	bufx := make([]float32, grid.nx)
	bufy := make([]float32, grid.ny)
	bufz := make([]float32, grid.nz)
	x := grid.GetVar("x")
	y := grid.GetVar("y")
	z := grid.GetVar("z")

	for i := 0; i < grid.nx; i++ {
		bufx[i] = float32(x[i*grid.di+0*grid.dj+0*grid.dk])
	}
	for j := 0; j < grid.ny; j++ {
		bufy[j] = float32(y[0*grid.di+j*grid.dj+0*grid.dk])
	}
	for k := 0; k < grid.nz; k++ {
		bufz[k] = float32(z[0*grid.di+0*grid.dj+k*grid.dk])
	}

	// fill buffer for boundary points
	ijk := 0
	for k := 0; k < grid.ny; k++ {
		for j := 0; j < grid.ny; j++ {
			for i := 0; i < grid.nx; i++ {
				fmt.Println(ijk, i*grid.di+j*grid.dj+k*grid.dk)
				buf[ijk] = int32(grid.boundary[i*grid.di+j*grid.dj+k*grid.dk])
				ijk++
			}
		}
	}

	r := f.Writer("boundary", d2, d1[1:grid.dim+1])
	r.Write(buf)
	r = f.Writer("x", []int{0}, []int{grid.nx})
	r.Write(bufx)
	if grid.dim > 1 {
		r = f.Writer("y", []int{0}, []int{grid.ny})
		r.Write(bufy)
	}
	if grid.dim > 2 {
		r = f.Writer("z", []int{0}, []int{grid.nz})
		r.Write(bufy)
	}

}

func (grid *Grid) Output_netcdf(it int) {

	fmt.Println("Write netcdf file:", outfile_netcdf, it)

	ff, _ := os.OpenFile(outfile_netcdf, os.O_RDWR, 660)
	f, _ := cdf.Open(ff)

	for _, vname := range grid.outvars {

		r := f.Writer(vname, []int{it, 0}, []int{it, grid.nx})
		switch grid.dim {
		case 2:
			r = f.Writer(vname, []int{it, 0, 0}, []int{it, grid.ny, grid.nx})
		case 3:
			r = f.Writer(vname, []int{it, 0, 0, 0}, []int{it, grid.nz, grid.ny, grid.nx})
		}

		buf := make([]float32, grid.nx*grid.ny*grid.nz) // a []T of the right T for the variable.

		v := grid.GetVar(vname)
		ijk := 0
		for k := 0; k < grid.nz; k++ {
			for j := 0; j < grid.ny; j++ {
				for i := 0; i < grid.nx; i++ {
					buf[ijk] = float32(v[i*grid.di+j*grid.dj+k*grid.dk])
					ijk++
				}
			}
		}

		fmt.Println(buf)

		r.Write(buf)
	}
	//os.Exit(0)
}
