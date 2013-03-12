package main

import (
	"fmt"
	"math"
	"os"

	. "github.com/marcusthierfelder/evolve"
)

var (
	outfile = os.Getenv("GOPATH") + "/src/out_adv.dat"
)

func rhs(grid *Grid, r, evl VarList) {
	f := evl.GetVar(0)
	df := r.GetVar(0)
	nx, dx := grid.GetSize()

	// simple symmetric stencil
	for i := 1; i < nx[0]-1; i++ {
		df[i] = -0.5 / dx[0] * (f[i+1] - f[i-1])

		//df[i] = -1./dx[0] * (f[i]-f[i-1])
	}

	// continius bc
	df[0], df[nx[0]-1] = df[nx[0]-2], df[1]
}

func main() {
	fmt.Println("Start an advection:")

	/* tell the library the domain and which variables do exist */
	grid := CreateGrid([]int{51}, []float64{0.}, []float64{1. / 50.})
	grid.AddVars([]string{"f"})

	/* create a variable list which is evolved */
	vl := grid.CreateVarlist()
	vl.AddVars([]string{"f"})

	/* initialize integrator */
	grid.TimeInt_init(vl, "rk4", interface{}(rhs))

	/* initial data */
	x := grid.GetVar("x")
	f := grid.GetVar("f")
	for i := 0; i < 51; i++ {
		f[i] = math.Exp(-50. * (x[i] - 0.5) * (x[i] - 0.5))
	}

	/* output file */
	fo, err := os.Create(outfile)
	if err != nil {
		panic(err)
	}
	defer fo.Close()
	grid.SetOutVars([]string{"f"})


	/* time step */
	dt := 0.01
	ot := 10
	grid.Output_netcdf_init(1000/ot)
	for t := 0; t < 1000; t++ {
		// output timestep
		if t%ot == 0 {
			fmt.Println("step", t, grid.GetTime())

			line := fmt.Sprintf("Time = %2.6f\n", grid.GetTime())
			fo.Write([]byte(line))
			for i, _ := range f {
				line = fmt.Sprintf("%2.6e %2.6e\n", x[i], f[i])
				fo.Write([]byte(line))
			}
			fo.Write([]byte("\n"))
			grid.Output_netcdf(t/ot)
		}
		// integrade one timestep
		grid.TimeInt(vl, &dt)
	}

}
