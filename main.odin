package main

import "core:fmt"
import mp_vec_sum "./more_precise_vec_sum"

main :: proc ( ) {

    fmt.printfln( "\nBegin a more precise vector sum in Odin ...\n\n" )

    mp_vec_sum.test_more_precise_vec_sum( )

    fmt.printfln( "\n\n... end of a more precise vector sum in Odin.\n" )
}
