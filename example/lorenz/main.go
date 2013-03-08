package main

import (
	"fmt"
	"os"

	. "github.com/marcusthierfelder/evolve"
)

var (
	outfile = os.Getenv("GOPATH") + "/src/out_lorenz.dat"

)

func rhs(grid *Grid, r, evl VarList) {
	x := evl.GetVar(0)
	y := evl.GetVar(1)
	z := evl.GetVar(2)

	dx := r.GetVar(0)
	dy := r.GetVar(1)
	dz := r.GetVar(2)

	sigma := 10.
	rho := 28.
	beta := 8./3.

	
	dx[0] = sigma*(y[0]-x[0])
	dy[0] = x[0]*(rho-z[0]) -y[0]
	dz[0] = x[0]*y[0] - beta*z[0]

}

func main() {
	fmt.Println("Start an advection:")

	/* tell the library the domain and which variables do exist */
	grid := CreateGrid([]int{}, []float64{}, []float64{})
	
	/* create a variable list which is evolved */
	vl := grid.CreateVarlist()
	vl.AddVars([]string{"x","y","z"})

	
	/* initialize integrator */
	grid.TimeInt_init(vl, "rk45", interface{}(rhs))

	/* initial data */
	x := grid.GetVar("x")
	y := grid.GetVar("y")
	z := grid.GetVar("z")
	x[0] = 1
	y[0] = 1
	z[0] = 1

	/* output file */
	fo, err := os.Create(outfile)
	if err != nil {
		panic(err)
	}
	defer fo.Close()

	/* time step */
	dt := 0.01
	ot := 10
	for t := 0; t < 10000; t++ {
		// output timestep
		if t%ot == 0 {
			fmt.Println("step", t, grid.GetTime())

			
			line := fmt.Sprintf("%2.6e %2.6e %2.6e\n", x[0], y[0],z[0])
			fo.Write([]byte(line))
			
		}
		// integrade one timestep
		grid.TimeInt(vl, &dt)
	}

}
