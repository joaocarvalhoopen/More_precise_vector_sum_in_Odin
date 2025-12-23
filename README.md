# More precise vector sum in Odin
Just a comparison of the some different high precision methods of making sums.

## Description
After the post of the last repository, I was currious to see what other interesting things would Prof. Cleve Moler had in is Blog, and I set my eyes on another article called SuperSum, In Defense of Floating Point Arithmetic by Cleve Moler. <br>
It's a very breaf article but one that set my curiosity on fire! <br>
I have read the article and then went to Wikipedia fantastic article on Kahan summation algorithm, or Compensated Sums, with many methods more sofisticated and more precise and faster the ones on Prof. Cleve Moler article.  And I follow that Wikipedia article links into the Pairwise Sums this last one is the one used by Pythons NumPy and Julia. I was hooked, so I looked into the core:math.odin file in the Odin code base for the sum( ) procedure and it is a simple sum function, then I converted some of the wikipedia examples and pseudocode of those 2 web pages into Odin, and implemented also an exact version with big.Int and big.Rat ( Rational ). Then I measured the time , so I have the exact value and the elapsed time of each execution. <br>
See the references below, and the results for an example for Odin below.

## Results

This are values for the sum of 1 Mega elements ( 1024 x 1024 ).


```
./mp_vec_sum_opti.exe

Begin a more precise vector sum in Odin ...


normal_sum                  = 55027917630045.828  Hex = 0h42C90617BDD52EEA  duration = 0.460 ms
                                       ^_____Big Error in value in Odin sum() !__________^

sum_long                    = 55027917629395.984  Hex = 0h42C90617BDD3E9FE  duration = 0.129 ms
                                         ^_sum_long() faster more precis then Odin sum()_^

kahan_sum                   = 55027917629453.828  Hex = 0h42C90617BDD406EA  duration = 1.646 ms
kahan_babushka_neumaier_sum = 55027917629453.828  Hex = 0h42C90617BDD406EA  duration = 0.790 ms
kahan_babushka_klein_sum    = 55027917629453.828  Hex = 0h42C90617BDD406EA  duration = 1.318 ms

shift_reduce_sum            = 55027917629453.828  Hex = 0h42C90617BDD406EA  duration = 0.517 ms

shift_reduce_sum_unrolled   = 55027917629453.828  Hex = 0h42C90617BDD406EA  duration = 0.886 ms

shift_reduce_sum_unrolled_2 = 55027917629453.828  Hex = 0h42C90617BDD406EA  duration = 0.520 ms
                                               Best Compromise!__________________________^

sum( f64 ) = 55027917629453.828  exact? = false
big_rat_rational_sum        = 55027917629453.828  Hex = 0h42C90617BDD406EA  duration = 1559.617 ms
big_sum exact_str           = 55027917629453827809570654574


... end of a more precise vector sum in Odin.

```

## References

1. SuperSum, In Defense of Floating Point Arithmetic <br>
   by Cleve Moler <br>
   [https://blogs.mathworks.com/cleve/2024/06/27/supersum-in-defense-of-floating-point-arithmetic/](https://blogs.mathworks.com/cleve/2024/06/27/supersum-in-defense-of-floating-point-arithmetic/)

2. Wikipedia - Kahan summation algorithm <br>
   [https://en.wikipedia.org/wiki/Kahan_summation_algorithm](https://en.wikipedia.org/wiki/Kahan_summation_algorithm)

3. Wikipedia - Pairwise summation <br>
   [https://en.wikipedia.org/wiki/Pairwise_summation](https://en.wikipedia.org/wiki/Pairwise_summation)

## How to compile and run

``` bash
make opti
make run_opti
```

## License
MIT Open Source License

## Have fun
Best regards, <br>
Joao Carvalho
