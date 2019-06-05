using MATLAB

x = range(-10.0, stop=10.0, length=500)
mat"plot($x, sin($x))"  # evaluate a MATLAB function

y = range(2.0, stop=3.0, length=500)
mat"""
    $u = $x + $y
	$v = $x - $y
"""
@show u v
