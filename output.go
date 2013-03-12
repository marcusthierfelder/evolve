package evolve

import (
	"os"
)

var (
	outpath = os.Getenv("GOPATH") + "/src"
	outfile = "out"
	outvars = []string{}
)

func (grid *Grid) SetOutVars(vars []string) {
	outvars = vars
}

func (grid *Grid) SetOutDir(dir, name string) {
	outpath = dir
	outfile = name
}
