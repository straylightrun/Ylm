// http://www.physics.drexel.edu/~tim/open/hydrofin/hyd.pdf

var a = 1; // Bohr Radius

// The coefficients used for computing in Legendre Polynomials
// http://arxiv.org/pdf/1410.1748.pdf
var A = []; // Coefficients Alm
var B = []; // Coefficients Blm
var P = []; // Plm
var Y = []; // Ylm
var LL = 64; // Max value for l

// Used sweet.js to compile this code to replace this macro
// PT(l,m) for all cases where you need to index A, B, P, 
/*
macro PT {
    rule {
        ($l:expr, $m:expr)
    } => {
         (( $m ) +(( $l ) *(( $l ) +1) ) /2)
    }
}

macro YR {
    rule {
        ($l:expr, $m:expr)
    } => {
         (( $m ) +( $l ) +(( $l ) *( $l ) ) )
    }
}
*/

function computeAB () {

	for(var l = 2; l < LL; l++) {
		var ls = l * l;
		var lm1s = (l - 1)*(l - 1);
		var ms = 0.0;

		for(var m = 0; m < l - 1; m++) {
			ms = m * m;
			A[m + l * (l + 1) / 2] = Math.sqrt((4 * ls - 1) / (ls - ms));
			B[m + l * (l + 1) / 2] = -Math.sqrt((lm1s - ms) / (4 * lm1s - 1));
		}
	}
}

function computeP(Lmax, x) {
	// Compute all Plm for given x
    var sintheta = Math.sqrt(1 - x * x);
    var temp = 0.3989422804014327; // = sqrt (0.5/ pi )
    P[0 + 0 * (0 + 1) / 2] = temp;
    if (Lmax > 0) {
        var sqrt3 = 1.7320508075688772;
        P[0 + 1 * (1 + 1) / 2] = x * sqrt3 * temp;
        sqrt3div2 = -1.224744871391589;
        temp = sqrt3div2 * sintheta * temp;
        P[1 + 1 * (1 + 1) / 2] = temp;
        for (var l = 2; l <= Lmax; l++) {
            for (var m = 0; m < l - 1; m++) {
                P[m + l * (l + 1) / 2] = A[m + l * (l + 1) / 2] * (x * P[m + (l - 1) * (l - 1 + 1) / 2] + B[m + l * (l + 1) / 2] * P[m + (l - 2) * (l - 2 + 1) / 2]);
            }
            P[l - 1 + l * (l + 1) / 2] = x * Math.sqrt(2 * (l - 1) + 3) * temp;
            temp = -Math.sqrt(1 + 0.5 / l) * sintheta * temp;
            P[l + l * (l + 1) / 2] = temp;
        }
    }
}

function computeY(Lmax, phi) {
    // Compute all Ylm for given phi
    var sqrt2 = 1.4142135623730951;
    var c1 = 1;
    var c2 = Math.cos(phi);
    var s1 = 0;
    var s2 = -Math.sin(phi);
    var tc = 2 * c2;
    for (var l = 0; l <= Lmax; l++) {
        Y[0 + l + l * l] = P[0 + l * (l + 1) / 2] * 0.5 * sqrt2;
    }
    for (var m = 1; m <= Lmax; m++) {
        var s = tc * s1 - s2;
        var c = tc * c1 - c2;
        s2 = s1;
        s1 = s;
        c2 = c1;
        c1 = c;
        for (var l = m; l <= Lmax; l++) {
            Y[-m + l + l * l] = P[m + l * (l + 1) / 2] * s;
            Y[m + l + l * l] = P[m + l * (l + 1) / 2] * c;
        }
    }
}

function L(l, n, r) {
	// Associated Laguerre polynomial
	// http://mathworld.wolfram.com/AssociatedLaguerrePolynomial.html

	var tmp = 0.0;
	for(var i = 0; i <= n; i++) {
		tmp += (math.factorial(n)/math.factorial(i)) * math.combinations(l + n, n - i) * math.pow(-r, i);
	}

	return (1/math.factorial(n)) * tmp;
}

function R(n, l, r) {
	// Radial Equation
	return math.sqrt(math.factorial(n - l - 1)/(2*n*(math.factorial(n+1)))) * math.exp(-4/(n*a)) * math.pow(((2*r)/(n*a)), l) * L(2*l + 1, n - l - 1, (2*r)/(n*a) );
}


function Y(l, m, theta, phi) {
	// Spherical Harmonics


}