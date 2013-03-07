package main

import (
    "fmt"
    . "github.com/marcusthierfelder/evolve"
//"reflect"
)




func (grid *Grid) rhs(r, evl VarList) {
    
}


func main() {
    fmt.Println("")


    grid := CreateGrid([]int{11}, []float64{0.}, []float64{1.})
    grid.AddVars([]string{"x","f"})


    x := grid.GetVar("x")

    f := grid.GetVar("x")

    for i:=0; i<11; i++ {
        x[i] = 1.        
    }
    f[0] = 1.








}