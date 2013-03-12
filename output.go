package evolve

import (
	"os"
)

var (
	outpath = os.Getenv("GOPATH") + "/src"
	outfile = "out"
	outvars = []string{}
)

// which variables have to be written to disk
func (grid *Grid) SetOutVars(vars []string) {
	outvars = vars
}

// set the path and the name of the output
func (grid *Grid) SetOutDir(dir, name string) {
	outpath = dir
	outfile = name
}
