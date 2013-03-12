package evolve

import (
	"fmt"
	"io"
	"log"
	"os"
)

func (grid *Grid) Output_vtk(it int) {

	/* init buffer */
	flag := os.O_CREATE | os.O_TRUNC | os.O_RDWR

	file := fmt.Sprintf("%s/%s_%04d.vtk", outpath, outfile, it)

	f, err := os.OpenFile(file, flag, 0666)
	if err != nil {
		log.Fatal(err)
		return
	}

	fmt.Println("Write vtk file: ", file, grid.time)
	io.WriteString(f, fmt.Sprintf("# vtk DataFile Version 2.0\n"))
	io.WriteString(f, fmt.Sprintf("variable f, level 0, time %16.9e\n", grid.time))
	io.WriteString(f, fmt.Sprintf("ASCII\n"))
	io.WriteString(f, fmt.Sprintf("DATASET STRUCTURED_POINTS\n"))
	io.WriteString(f, fmt.Sprintf("DIMENSIONS %d %d %d\n", grid.nx, grid.ny, grid.nz))
	io.WriteString(f, fmt.Sprintf("ORIGIN  %16.9e %16.9e %16.9e\n", grid.x0, grid.y0, grid.z0))
	io.WriteString(f, fmt.Sprintf("SPACING %16.9e %16.9e %16.9e\n", grid.dx, grid.dy, grid.dz))
	io.WriteString(f, fmt.Sprintf("\n"))
	io.WriteString(f, fmt.Sprintf("POINT_DATA %d\n", grid.nx*grid.ny*grid.nz))

	for _, v := range outvars {
		d := grid.GetVar(v)

		io.WriteString(f, fmt.Sprintf("SCALARS %s float\n", v))
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
	}

	f.Close()
}
