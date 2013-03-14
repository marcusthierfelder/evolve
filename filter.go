package evolve

import (
	_ "fmt"
	"math"
	_ "os"
	_ "reflect"
)

func (grid *Grid) damping(r, uc VarList) {

	NX, DI, DX := grid.GetSize()
	nx, ny := NX[0], NX[1]
	di, dj := DI[0], DI[1]
	dx, dy := DX[0], DX[1]

	dissfactor := 0.2
	oo2dx := 1. / (2. * dx)
	oo2dy := 1. / (2. * dy)

	for iv, _ := range r.field {
		v := uc.GetVar(iv)
		r := r.GetVar(iv)

		switch grid.dim {
		case 2:
			for j := 0; j < ny; j++ {
				for i := 0; i < nx; i++ {
					ij := j*dj + i*di
					// continue

					d4 := 0.
					if i > 1 && i < nx-2 && j > 1 && j < ny-2 {
						d4 = -2. * dissfactor * (oo2dx*(6.*v[ij]+v[-2*di+ij]-4.*(v[-di+ij]+v[di+ij])+v[2*di+ij]) +
							oo2dy*(6.*v[ij]+v[-2*dj+ij]-4.*(v[-dj+ij]+v[dj+ij])+v[2*dj+ij]))

						if math.IsNaN(d4) {
							d4 = 0.
						}
					}

					d2 := 0.
					if i > 0 && i < nx-1 && j > 0 && j < ny-1 {
						d2 = -2. * dissfactor * (oo2dx*(-2.*v[ij]+1.*(v[-di+ij]+v[di+ij])) +
							oo2dy*(-2.*v[ij]+1.*(v[-dj+ij]+v[dj+ij])))

						if math.IsNaN(d2) {
							d2 = 0.
						}
					}

					//fmt.Println(d4)

					r[ij] -= 0.*d2 + d4

				}
			}

		default:
			panic("damping for this dimention is not implemented")

		}

	}

}
