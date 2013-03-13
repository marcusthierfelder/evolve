package evolve

import (
	"os"
)

var (
	outpath = os.Getenv("GOPATH") + "/src"
	outfile = "out"
)

// which variables have to be written to disk
func (grid *Grid) SetOutVars(vars []string) {
	grid.outvars = vars
}

// set the path and the name of the output
func (grid *Grid) SetOutDir(dir, name string) {
	outpath = dir
	outfile = name
}
