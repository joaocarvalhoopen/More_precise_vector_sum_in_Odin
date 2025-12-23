package more_precise_vec_sum

import "core:fmt"
import "core:math"
import "base:intrinsics"
import "core:math/big"
import "core:os"
import "core:strings"
import "core:time"

// from core/math.odin line 1629
@(require_results)
sum :: proc "contextless" (x: $T/[]$E) -> (res: E)
	where intrinsics.type_is_numeric(E) {
	for i in x {
		res += i
	}
	return
}

// Original Kahan Summation Algorithm - Compensated sum.
// https://en.wikipedia.org/wiki/Kahan_summation_algorithm
@(require_results)
kahan_sum :: proc "contextless" ( x : $T / [ ]$E ) -> E
	where intrinsics.type_is_numeric( E ) {

    // Prepare the accumulator.
    sum : E = 0.0

    // A running compensation for lost low-order bits.
    c   : E = 0.0

    for val in x {

        // c is zero the first time around.
        y := val - c
        // Alas, sum is big, y small, so low-order digits of y are lost.
        t := sum + y
        // ( t - sum ) cancels the high-order part of y;
        // subtracting y recovers negative ( low part of y )
        c = ( t - sum ) - y
        // Algebraically, c should always be zero. Beware
        // overly-aggressive optimizing compilers!
        sum = t
        // Next time around, the lost low part will be added to y in a fresh attempt.
    }

    return sum
}

// Improved Neumaier.
// https://en.wikipedia.org/wiki/Kahan_summation_algorithm
@(require_results)
kahan_babushka_neumaier_sum :: proc "contextless" ( x : $T / [ ]$E ) -> E
	where intrinsics.type_is_numeric( E ) {

    sum : E = 0.0
    c   : E = 0.0    // A running compensation for lost low-order bits.

    for val in x {

        t := sum + val
        if abs( sum ) >= abs( val ) {

            c += ( sum - t ) + val   // If sum is bigger, low-order digits of val are lost.
        } else {

            c += ( val - t ) + sum   // Else low-order digits of sum are lost.
        }
        sum = t
    }

    return sum + c   // Correction only applied once in the very end.
}

// Better accuracy then Neumaier.
// https://en.wikipedia.org/wiki/Kahan_summation_algorithm
@(require_results)
kahan_babushka_klein_sum :: proc "contextless" ( x : $T / [ ]$E ) -> E
	where intrinsics.type_is_numeric( E ) {

    sum : E = 0.0
    cs  : E = 0.0
    ccs : E = 0.0

    for val in x {

        c  : E
        cc : E
        t  : E = sum + val

        if abs( sum ) >= abs( val ) {

            c = ( sum - t ) + val
        } else {

            c = ( val - t ) + sum
        }

        sum = t
        t = cs + c

        if abs( cs ) >= abs(c ) {

            cc = ( cs - t ) + c
        } else {

            cc = ( c - t ) + cs
        }

        cs = t
        ccs = ccs + cc
	}

    return sum + ( cs + ccs )
}

// Pairwise Sumation
// https://en.wikipedia.org/wiki/Pairwise_summation
@(require_results)
shift_reduce_sum :: proc "contextless" ( x : $T / [ ]$E ) -> E
    where intrinsics.type_is_numeric( E ) {


    n : int = len( x )

    stack : [ 64 ]E
    v     : E
    p     : uint = 0

    for i in uint( 0 ) ..< uint( n ) {

        v = x[ i ]                                   // shift
        // for b : int = 1; i & b; b <<= 1, −−p      // reduce
        for b : uint = 1; i & b != 0; {              // reduce

            v += stack[ p - 1 ]
            b = b << 1
            p -= 1
        }

        // stack[ p++ ] = v
        stack[ p ] = v
        p += 1
    }

    sum : E = 0.0
    for p != 0 {

        p   -= 1
        sum += stack[ p ]
        // sum += stack[−−p];
    }

    return sum
}


/*

double shift_reduce_sum(double ∗x, size_t n) {
  double stack[64], v;
  size_t p = 0;
  for (size_t i = 0; i < n; ++i) {
    v = x[i];                               // shift
    for (size_t b = 1; i & b; b <<= 1, −−p) // reduce
      v += stack[p−1];
    stack[p++] = v;
  }
  double sum = 0.0;
  while (p)
    sum += stack[−−p];
  return sum;
}

*/





/*

// Big Int
@(require_results)
big_int_sum :: proc ( x : [ ]f64 ) -> f64 {

    a, b, sum: big.Int

    // Allocate backing storage for multiple Ints at once
	if err := big.int_init_multi( & a, & b, & sum ); err != .None {

	    fmt.println( "init error:", err )
		os.exit( -1 )
	}
	defer big.int_destroy( & a, & b, & sum )


	for elem in x {


	    _ = big.int_set_from_integer( & a, 123456789 )
		// _ = big.int_set_from_integer( & b, 987654321 )
		_ = big.int_add( & sum, & a, & sum )


		/*
		_ = big.int_set_from_integer( & a, 123456789 )
		_ = big.int_set_from_integer( & b, 987654321 )
		_ = big.int_add( & sum, & a, & b )
		*/
	}



	// int_to_string allocates; caller must free
	s, _ := big.int_to_string( & sum )
	defer delete( s )
	fmt.println( "sum =", s )

	return f64(  )
}

*/


// Big Rational
@(require_results)
big_rat_rational_sum :: proc ( x : [ ]f64 ) -> ( f64, string ) {

    sum : big.Rat
    tmp : big.Rat

	// Inicializa os Int internos ( sum.a, sum.b, tmp.a, tmp.b )
	if err := big.int_init_multi( & sum.a, & sum.b, & tmp.a, & tmp.b ); err != .None {

	    fmt.println( "init err:", err )
		os.exit( -1 )
	}
	defer big.int_destroy( & sum.a, & sum.b, & tmp.a, & tmp.b )

	// sum = 0
	_ = big.rat_set_u64( & sum, 0 )

	for v in x {

		// tmp = v (como racional)
		if err := big.rat_set_f64( & tmp, v ); err != .None {

			fmt.println( "rat_set_f64 err:", err )
			os.exit( -1 )
		}

		// sum = sum + tmp
		if err := big.rat_add_rat( & sum, & sum, & tmp ); err != .None {

			fmt.println("rat_add_rat err:", err)
			os.exit( -1 )
		}
	}

	// Voltar para f64
	f, exact, err := big.rat_to_f64( & sum )
	if err != .None {

		fmt.println( "rat_to_f64 err:", err )
		os.exit( -1 )
	}

	fmt.printfln( "sum( f64 ) = %17f  exact? = %v" , f , exact )


	// NOTE:
	//    In Big Decimal ( big.Rat ) and Big.Int.
	//    632 decimal digits holds every IEEE-754 f64
    //    632 = math.ceil( math.log10( real_max ) - math.log10( math.f64_epsilon * real_min ) )


    // Obtain the result has a single integer to have more precision then 16 decimal places.
	result : big.Int
	mul_factor : big.Int

	if err := big.int_init_multi( & result, & mul_factor ); err != .None {

	    fmt.println( "init result err:", err )
		os.exit( -1 )
	}
	defer big.int_destroy( & result )

	_ = big.int_set_from_integer( & mul_factor, 1_000_000_000_000_000 )

	_ = big.int_mul( & sum.a, & sum.a, & mul_factor )


	_ = big.int_div( & result, & sum.a, & sum.b )

	result_str, err_2 := big.int_to_string( & result )

	return f, result_str
}

/*

rat_to_string :: proc( r : ^big.Rat, allocator := context.allocator ) -> ( s : string, ok : bool ) {

    // numerator
    num, err1 := big.int_itoa_string( & r.a, allocator = allocator ) // radix=10 default
    if err1 != nil {

        return "", false
    }
    defer delete( num )

    // denominator
    den, err2 := big.int_itoa_string( & r.b, allocator = allocator )
    if err2 != nil {

        return "", false
    }
    defer delete( den )

    // nice formatting: if denominator is 1, just return numerator
    if den == "1" {

        // need a clone because num will be deleted on defer
        out, aerr := strings.concatenate( [ ? ]string{ num }[ : ], allocator = allocator )
        return out, aerr == nil
    }

    out, aerr := strings.concatenate( [ ? ]string{ num, "/", den }[ : ], allocator = allocator )
    return out, aerr == nil
}

*/


test_more_precise_vec_sum :: proc ( ) {

    VEC_LEN :: 1 * 1024 * 1024
    vec_a : [ ]f64 = make( [ ]f64, VEC_LEN )
    defer delete( vec_a )

    // Fill the vec_a .
    val : f64 = 0.123456789012345678
    for & elem, i in vec_a {

        elem = val
        val += math.F64_EPSILON * 1.111_111_111_111_111_111

        if i % 1000 == 0 {

            val += 1e5
        }
    }


   	start := time.now()

    res_normal_sum : f64 = sum( vec_a )

    elapsed := time.since( start )
    ms      := time.duration_milliseconds( elapsed )

    fmt.printfln( "normal_sum                  = %17f  Hex = %H  duration = %.3f ms",
                  res_normal_sum, res_normal_sum, ms )
    fmt.printfln( "                                       ^_____Big Error in value in Odin sum() !\n" )


    start = time.now()

    res_kahan_sum : f64 = kahan_sum( vec_a )

    elapsed = time.since( start )
    ms      = time.duration_milliseconds( elapsed )


    fmt.printfln( "kahan_sum                   = %17f  Hex = %H  duration = %.3f ms",
                  res_kahan_sum, res_kahan_sum, ms )


    start = time.now()

    res_kahan_babushka_neumaier_sum : f64 = kahan_babushka_neumaier_sum( vec_a )

    elapsed = time.since( start )
    ms      = time.duration_milliseconds( elapsed )

    fmt.printfln( "kahan_babushka_neumaier_sum = %17f  Hex = %H  duration = %.3f ms",
                  res_kahan_babushka_neumaier_sum, res_kahan_babushka_neumaier_sum, ms )


    start = time.now()

    res_kahan_babushka_klein_sum : f64 = kahan_babushka_klein_sum( vec_a )

    elapsed = time.since( start )
    ms      = time.duration_milliseconds( elapsed )

    fmt.printfln( "kahan_babushka_klein_sum    = %17f  Hex = %H  duration = %.3f ms\n",
                  res_kahan_babushka_klein_sum, res_kahan_babushka_klein_sum, ms )


    start = time.now()

    res_shift_reduce_sum : f64 = shift_reduce_sum( vec_a )

    elapsed = time.since( start )
    ms      = time.duration_milliseconds( elapsed )

    fmt.printfln( "shift_reduce_sum            = %17f  Hex = %H  duration = %.3f ms\n",
                  res_shift_reduce_sum, res_shift_reduce_sum, ms )

    start = time.now()

    res_big_rat_rational_sum, result_exact_str := big_rat_rational_sum( vec_a )

    elapsed = time.since( start )
    ms      = time.duration_milliseconds( elapsed )

    fmt.printfln( "big_rat_rational_sum        = %17f  Hex = %H  duration = %.3f ms",
                   res_big_rat_rational_sum, res_big_rat_rational_sum, ms )
    fmt.printfln( "big_sum exact_str           = %s", result_exact_str )
}
