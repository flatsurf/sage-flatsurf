r"""
Utilities to deal with power series and Laurent series defined on surfaces.
"""
######################################################################
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2022-2023 Julian Rüth
#
#  sage-flatsurf is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 2 of the License, or
#  (at your option) any later version.
#
#  sage-flatsurf is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with sage-flatsurf. If not, see <https://www.gnu.org/licenses/>.
######################################################################
from sage.rings.ring import CommutativeRing
from sage.structure.element import CommutativeRingElement
from sage.structure.unique_representation import UniqueRepresentation
from sage.misc.cachefunc import cached_method


class PowerSeriesCoefficientExpression(CommutativeRingElement):
    r"""
    An expression in the (symbolic) coefficients of a multivariate power
    series.

    Consider a multivariate power series, for example, in two variables

    .. MATH::

        \sum a_i b_j x^i y^j

    This element is an expression in the coefficients of such a series, for
    example

    .. MATH::

        a_0^2 + 2 a_1 b_1 + b_1^2

    In principle, this is just a multivariate polynomial in all the `a_i` and
    `b_j`. However, for our use case, see :module:`harmonic_differentials`,
    these expressions are of low degree (at most 2) but with lots of variables.
    Multivariate polynomial rings in SageMath are bad at handling such
    extremely sparse scenarios since some operations are implemented linearly
    in the number of generators of the ring.

    Therefore, we roll our own "sparse multivariate polynomial ring" here that
    is asymptotically fast for the operations we care about (and slow for other
    operations such as multiplication of expressions.)

    EXAMPLES::

        sage: from flatsurf.geometry.power_series import PowerSeriesCoefficientExpressionRing  # random output due to deprecation warnings from cppyy
        sage: R = PowerSeriesCoefficientExpressionRing(QQ, gens=("a_?", "b_?"))
        sage: a0 = R.gen(("a_?", 0))
        sage: b0 = R.gen(("b_?", 0))
        sage: a1 = R.gen(("a_?", 1))
        sage: b1 = R.gen(("b_?", 1))

        sage: a0^2 + b0^2 + 2 * a1 * b1
        a_0^2 + b_0^2 + 2*a_1*b_1

    """

    def __init__(self, parent, coefficients):
        super().__init__(parent)

        # Zero coefficients must be removed. Otherwise, degree computations break.
        assert all(v for v in coefficients.values())

        self._coefficients = coefficients

    def _richcmp_(self, other, op):
        r"""
        Compare this expression to ``other``.

        EXAMPLES::

            sage: from flatsurf.geometry.power_series import PowerSeriesCoefficientExpressionRing
            sage: R = PowerSeriesCoefficientExpressionRing(QQ, gens=("a_?", "b_?"))
            sage: a0 = R.gen(("a_?", 0))
            sage: b0 = R.gen(("b_?", 0))

            sage: a0 == b0
            False

            sage: a0^2 == a0
            False

            sage: a0 == a0
            True

        """
        from sage.structure.richcmp import op_EQ, op_NE

        if op == op_NE:
            return not (self == other)

        if op == op_EQ:
            return self._coefficients == other._coefficients

        raise NotImplementedError

    def _repr_(self):
        r"""
        Return a printable representation of this expression.

        EXAMPLES::

            sage: from flatsurf.geometry.power_series import PowerSeriesCoefficientExpressionRing
            sage: R = PowerSeriesCoefficientExpressionRing(QQ, gens=("Re(a0,?)", "Im(a0,?)"))
            sage: re = R.gen(("Re(a0,?)", 0))
            sage: im = R.gen(("Im(a0,?)", 0))

            sage: re
            Re(a0,0)
            sage: im
            Im(a0,0)
            sage: re + im
            Re(a0,0) + Im(a0,0)
            sage: re + im + 1
            Re(a0,0) + Im(a0,0) + 1

        """
        def variable_name(variable):
            gen, degree = variable.describe()

            return gen.replace("?", str(degree))

        variables = list(self.variables())

        def key(variable):
            return next(iter(next(iter(variable._coefficients.keys()))))

        variables.sort(key=key)

        variable_names = tuple(variable_name(variable) for variable in variables)

        def encode_variable_name(variable):
            return variable.replace("(", "__open__").replace(")", "__close__").replace(",", "__comma__")

        variable_names = [encode_variable_name(name) for name in variable_names]

        def decode_variable_name(variable):
            return variable.replace("__comma__", ",").replace("__close__", ")").replace("__open__", "(")

        from sage.all import PolynomialRing
        R = PolynomialRing(self.base_ring(), variable_names)

        def monomial(gens):
            monomial = R.one()
            for gen in gens:
                gen = self.parent().gen(gen)
                monomial *= R(variable_names[variables.index(gen)])
            return monomial

        f = sum(coefficient * monomial(gens) for (gens, coefficient) in self._coefficients.items())

        return decode_variable_name(repr(f))

    def degree(self, gen):
        r"""
        Return the total degree of this expression in the variable ``gen``.

        EXAMPLES::

            sage: from flatsurf.geometry.power_series import PowerSeriesCoefficientExpressionRing
            sage: R = PowerSeriesCoefficientExpressionRing(QQ, gens=("Re(a0,?)", "Im(a0,?)"))
            sage: a = R.gen(("Re(a0,?)", 0))
            sage: b = R.gen(("Im(a0,?)", 0))

            sage: a.degree(a)
            1
            sage: (a + b).degree(a)
            1
            sage: (a * b + a).degree(a)
            1
            sage: R.one().degree(a)
            0
            sage: R.zero().degree(a)
            -1

        """
        if not gen.is_variable():
            raise ValueError(f"gen must be a variable not {gen}")

        variable = next(iter(next(iter(gen._coefficients))))

        return max([monomial.count(variable) for monomial in self._coefficients], default=-1)

    def is_monomial(self):
        r"""
        Return whether this expression is a non-constant monomial without a
        leading coefficient.

        EXAMPLES::

            sage: from flatsurf.geometry.power_series import PowerSeriesCoefficientExpressionRing
            sage: R = PowerSeriesCoefficientExpressionRing(QQ, gens=("Re(a0,?)", "Im(a0,?)"))
            sage: a = R.gen(("Re(a0,?)", 0))
            sage: b = R.gen(("Im(a0,?)", 0))

            sage: a.is_monomial()
            True
            sage: (2*a).is_monomial()
            False
            sage: (a + b).is_monomial()
            False
            sage: R.one().is_monomial()
            False
            sage: R.zero().is_monomial()
            False
            sage: (a * a).is_monomial()
            True

        """
        if len(self._coefficients) != 1:
            return False

        ((key, value),) = self._coefficients.items()

        return bool(key) and value.is_one()

    def is_constant(self):
        r"""
        Return whether this expression is a constant from the base ring.

        EXAMPLES::

            sage: from flatsurf.geometry.power_series import PowerSeriesCoefficientExpressionRing
            sage: R = PowerSeriesCoefficientExpressionRing(QQ, gens=("Re(a0,?)", "Im(a0,?)"))
            sage: a = R.gen(("Re(a0,?)", 0))
            sage: b = R.gen(("Im(a0,?)", 0))

            sage: a.is_constant()
            False
            sage: (a + b).is_constant()
            False
            sage: R.one().is_constant()
            True
            sage: R.zero().is_constant()
            True
            sage: (a * a).is_constant()
            False

        """
        coefficients = len(self._coefficients)

        if coefficients == 0:
            return True

        if coefficients > 1:
            return False

        monomial = next(iter(self._coefficients.keys()))

        return not monomial

    def norm(self, p=2):
        r"""
        Return the p-norm of the coefficient vector.

        EXAMPLES::

            sage: from flatsurf.geometry.power_series import PowerSeriesCoefficientExpressionRing
            sage: R = PowerSeriesCoefficientExpressionRing(QQ, gens=("Re(a0,?)", "Im(a0,?)"))
            sage: a = R.gen(("Re(a0,?)", 0))
            sage: b = R.gen(("Im(a0,?)", 0))

            sage: x = a + 1

            sage: x.norm(1)
            2

            sage: x.norm(oo)
            1

        """
        from sage.all import vector
        return vector(self._coefficients.values()).norm(p)

    def _neg_(self):
        r"""
        EXAMPLES::

            sage: from flatsurf.geometry.power_series import PowerSeriesCoefficientExpressionRing
            sage: R = PowerSeriesCoefficientExpressionRing(QQ, gens=("Re(a0,?)", "Im(a0,?)"))
            sage: a = R.gen(("Re(a0,?)", 0))
            sage: b = R.gen(("Im(a0,?)", 0))

            sage: -a
            -Re(a0,0)
            sage: -(a + b)
            -Re(a0,0) - Im(a0,0)
            sage: -(a * a)
            -Re(a0,0)^2
            sage: -R.one()
            -1
            sage: -R.zero()
            0

        """
        parent = self.parent()
        return type(self)(parent, {key: -coefficient for (key, coefficient) in self._coefficients.items()})

    def _add_(self, other):
        r"""
        EXAMPLES::

            sage: from flatsurf.geometry.power_series import PowerSeriesCoefficientExpressionRing
            sage: R = PowerSeriesCoefficientExpressionRing(QQ, gens=("Re(a0,?)", "Im(a0,?)"))
            sage: a = R.gen(("Re(a0,?)", 0))
            sage: b = R.gen(("Im(a0,?)", 0))

            sage: a + 1
            Re(a0,0) + 1
            sage: a + (-a)
            0
            sage: a + b
            Re(a0,0) + Im(a0,0)
            sage: a * a + b * b
            Re(a0,0)^2 + Im(a0,0)^2

        """
        return self.parent().sum([self, other])

    def _sub_(self, other):
        r"""
        EXAMPLES::

            sage: from flatsurf.geometry.power_series import PowerSeriesCoefficientExpressionRing
            sage: R = PowerSeriesCoefficientExpressionRing(QQ, gens=("Re(a0,?)", "Im(a0,?)"))
            sage: a = R.gen(("Re(a0,?)", 0))
            sage: b = R.gen(("Im(a0,?)", 0))

            sage: a - 1
            Re(a0,0) - 1
            sage: a - a
            0
            sage: a * a - b * b
            Re(a0,0)^2 - Im(a0,0)^2

        """
        return self._add_(-other)

    def _mul_(self, other):
        r"""
        EXAMPLES::

            sage: from flatsurf.geometry.power_series import PowerSeriesCoefficientExpressionRing
            sage: R = PowerSeriesCoefficientExpressionRing(QQ, gens=("Re(a0,?)", "Im(a0,?)"))
            sage: a = R.gen(("Re(a0,?)", 0))
            sage: b = R.gen(("Im(a0,?)", 0))

            sage: a * a
            Re(a0,0)^2
            sage: a * b
            Re(a0,0)*Im(a0,0)
            sage: a * R.one()
            Re(a0,0)
            sage: a * R.zero()
            0
            sage: (a + b) * (a - b)
            Re(a0,0)^2 - Im(a0,0)^2

        """
        parent = self.parent()

        if other.is_zero() or self.is_zero():
            return parent.zero()

        if other.is_one():
            return self

        if self.is_one():
            return other

        coefficients = {}

        for self_monomial, self_coefficient in self._coefficients.items():
            assert self_coefficient
            for other_monomial, other_coefficient in other._coefficients.items():
                assert other_coefficient

                monomial = tuple(sorted(self_monomial + other_monomial))
                coefficient = self_coefficient * other_coefficient

                if monomial not in coefficients:
                    coefficients[monomial] = coefficient
                else:
                    coefficients[monomial] += coefficient
                    if not coefficients[monomial]:
                        del coefficients[monomial]

        return type(self)(self.parent(), coefficients)

    def _rmul_(self, right):
        r"""
        EXAMPLES::

            sage: from flatsurf.geometry.power_series import PowerSeriesCoefficientExpressionRing
            sage: R = PowerSeriesCoefficientExpressionRing(QQ, gens=("Re(a0,?)", "Im(a0,?)"))
            sage: a = R.gen(("Re(a0,?)", 0))
            sage: b = R.gen(("Im(a0,?)", 0))

            sage: a * 0
            0
            sage: a * 1
            Re(a0,0)
            sage: a * 2
            2*Re(a0,0)

        """
        return self._lmul_(right)

    def _lmul_(self, left):
        r"""
        EXAMPLES::

            sage: from flatsurf.geometry.power_series import PowerSeriesCoefficientExpressionRing
            sage: R = PowerSeriesCoefficientExpressionRing(QQ, gens=("Re(a0,?)", "Im(a0,?)"))
            sage: a = R.gen(("Re(a0,?)", 0))
            sage: b = R.gen(("Im(a0,?)", 0))

            sage: 0 * a
            0
            sage: 1 * a
            Re(a0,0)
            sage: 2 * a
            2*Re(a0,0)

        """
        return type(self)(self.parent(), {key: coefficient for (key, value) in self._coefficients.items() if (coefficient := left * value)})

    def constant_coefficient(self):
        r"""
        Return the constant coefficient of this expression.

        EXAMPLES::

            sage: from flatsurf.geometry.power_series import PowerSeriesCoefficientExpressionRing
            sage: R = PowerSeriesCoefficientExpressionRing(QQ, gens=("Re(a0,?)", "Im(a0,?)"))
            sage: a = R.gen(("Re(a0,?)", 0))
            sage: b = R.gen(("Im(a0,?)", 0))

            sage: a.constant_coefficient()
            0
            sage: (a + b).constant_coefficient()
            0
            sage: R.one().constant_coefficient()
            1
            sage: R.zero().constant_coefficient()
            0

        """
        return self._coefficients.get((), self.parent().base_ring().zero())

    def variables(self):
        r"""
        Return the variables that appear in this expression.

        EXAMPLES::

            sage: from flatsurf.geometry.power_series import PowerSeriesCoefficientExpressionRing
            sage: R = PowerSeriesCoefficientExpressionRing(QQ, gens=("Re(a0,?)", "Im(a0,?)"))
            sage: re = R.gen(("Re(a0,?)", 0))
            sage: im = R.gen(("Im(a0,?)", 0))

            sage: R.zero().variables()
            set()
            sage: R.one().variables()
            set()
            sage: (re^2 * im + 1).variables()
            {Im(a0,0), Re(a0,0)}
            sage: (re + 1).variables()
            {Re(a0,0)}

        """
        return set(self.parent().gen(gen) for monomial in self._coefficients.keys() for gen in monomial)

    def is_variable(self):
        r"""
        EXAMPLES::

            sage: from flatsurf.geometry.power_series import PowerSeriesCoefficientExpressionRing
            sage: R = PowerSeriesCoefficientExpressionRing(QQ, gens=("Re(a0,?)", "Im(a0,?)"))
            sage: re = R.gen(("Re(a0,?)", 0))
            sage: im = R.gen(("Im(a0,?)", 0))

            sage: re.is_variable()
            True
            sage: R.zero().is_variable()
            False
            sage: R.one().is_variable()
            False
            sage: (re + 1).is_variable()
            False
            sage: (re * im).is_variable()
            False

        """
        if not self.is_monomial():
            return False

        monomial = next(iter(self._coefficients.keys()))

        if len(monomial) != 1:
            return False

        return True

    def describe(self):
        r"""
        Return a tuple describing the nature of this variable.

        EXAMPLES::

            sage: from flatsurf.geometry.power_series import PowerSeriesCoefficientExpressionRing
            sage: R = PowerSeriesCoefficientExpressionRing(QQ, gens=("Re(a0,?)", "Im(a0,?)"))
            sage: re = R.gen(("Re(a0,?)", 0))
            sage: im = R.gen(("Im(a0,?)", 0))

            sage: re.describe()
            ('Re(a0,?)', 0)
            sage: im.describe()
            ('Im(a0,?)', 0)
            sage: (re + im).describe()
            Traceback (most recent call last):
            ...
            ValueError: element must be a variable

        """
        if not self.is_variable():
            raise ValueError("element must be a variable")

        variable = next(iter(next(iter(self._coefficients.keys()))))

        gens = self.parent()._gens

        gen = gens[variable % len(gens)]
        degree = variable // len(gens)

        return gen, degree

    def real(self):
        r"""
        Return the real part of this expression.

        EXAMPLES::

            sage: from flatsurf.geometry.power_series import PowerSeriesCoefficientExpressionRing
            sage: R = PowerSeriesCoefficientExpressionRing(CC, gens=("Re(a0,?)", "Im(a0,?)"))
            sage: re = R.gen(("Re(a0,?)", 0))
            sage: im = R.gen(("Im(a0,?)", 0))

            sage: x = (re + I*im)**2
            sage: x
            Re(a0,0)^2 + 2.00000000000000*I*Re(a0,0)*Im(a0,0) - Im(a0,0)^2
            sage: x.real()
            Re(a0,0)^2 - Im(a0,0)^2

        """
        return self.map_coefficients(lambda c: c.real(), self.parent().change_ring(self.parent().real_field()))

    def imag(self):
        r"""
        Return the imaginary part of this expression.

        EXAMPLES::

            sage: from flatsurf.geometry.power_series import PowerSeriesCoefficientExpressionRing
            sage: R = PowerSeriesCoefficientExpressionRing(CC, gens=("Re(a0,?)", "Im(a0,?)"))
            sage: re = R.gen(("Re(a0,?)", 0))
            sage: im = R.gen(("Im(a0,?)", 0))

            sage: x = (re + I*im)**2
            sage: x
            Re(a0,0)^2 + 2.00000000000000*I*Re(a0,0)*Im(a0,0) - Im(a0,0)^2
            sage: x.imag()
            2.00000000000000*Re(a0,0)*Im(a0,0)

        """
        return self.map_coefficients(lambda c: c.imag(), self.parent().change_ring(self.parent().real_field()))

    def __getitem__(self, gen):
        r"""
        Return the coefficient of the monomial ``gen``.

        EXAMPLES::

            sage: from flatsurf.geometry.power_series import PowerSeriesCoefficientExpressionRing
            sage: R = PowerSeriesCoefficientExpressionRing(CC, gens=("Re(a0,?)", "Im(a0,?)"))
            sage: re = R.gen(("Re(a0,?)", 0))
            sage: im = R.gen(("Im(a0,?)", 0))

            sage: re[re]
            1.00000000000000
            sage: re[im]
            0.000000000000000
            sage: (re + im)[re]
            1.00000000000000
            sage: (re * im)[re]
            0.000000000000000
            sage: (re * im)[re * im]
            1.00000000000000

        """
        if not gen.is_monomial():
            raise ValueError("gen must be a monomial")

        return self._coefficients.get(next(iter(gen._coefficients.keys())), self.parent().base_ring().zero())

    def __hash__(self):
        return hash(tuple(sorted(self._coefficients.items())))

    def total_degree(self):
        r"""
        Return the total degree of this expression in its symbolic variables.

        EXAMPLES::

            sage: from flatsurf.geometry.power_series import PowerSeriesCoefficientExpressionRing
            sage: R = PowerSeriesCoefficientExpressionRing(CC, gens=("Re(a0,?)", "Im(a0,?)"))
            sage: re = R.gen(("Re(a0,?)", 0))
            sage: im = R.gen(("Im(a0,?)", 0))

            sage: R.zero().total_degree()
            -1
            sage: R.one().total_degree()
            0
            sage: re.total_degree()
            1
            sage: (re * re + im).total_degree()
            2

        """
        degrees = [len(monomial) for monomial in self._coefficients]
        return max(degrees, default=-1)

    def derivative(self, gen):
        r"""
        Return the derivative of this expression with respect to the variable
        ``gen``.

        EXAMPLES::

            sage: from flatsurf.geometry.power_series import PowerSeriesCoefficientExpressionRing
            sage: R = PowerSeriesCoefficientExpressionRing(CC, gens=("Re(a0,?)", "Im(a0,?)"))
            sage: re = R.gen(("Re(a0,?)", 0))
            sage: im = R.gen(("Im(a0,?)", 0))

            sage: R.zero().derivative(re)
            0
            sage: R.one().derivative(re)
            0
            sage: re.derivative(re)
            1.00000000000000
            sage: re.derivative(im)
            0
            sage: x = (re + im) * (re - im)
            sage: x.derivative(re)
            2.00000000000000*Re(a0,0)
            sage: x.derivative(im)
            -2.00000000000000*Im(a0,0)

        """
        if not gen.is_variable():
            raise ValueError("gen must be a variable")

        gen = next(iter(gen._coefficients.keys()))[0]

        derivative = self.parent().zero()

        for monomial, coefficient in self._coefficients.items():
            assert coefficient

            exponent = monomial.count(gen)

            if not exponent:
                continue

            monomial = list(monomial)
            monomial.remove(gen)
            monomial = tuple(monomial)

            derivative += self.parent()({monomial: exponent * coefficient})

        return derivative

    def map_coefficients(self, f, ring=None):
        r"""
        Return the image of this expression by applying ``f`` to each non-zero
        coefficient of the expression.

        EXAMPLES::

            sage: from flatsurf.geometry.power_series import PowerSeriesCoefficientExpressionRing
            sage: R = PowerSeriesCoefficientExpressionRing(CC, gens=("Re(a0,?)", "Im(a0,?)"))
            sage: re = R.gen(("Re(a0,?)", 0))
            sage: im = R.gen(("Im(a0,?)", 0))

            sage: re.map_coefficients(lambda c: 2*c)
            2.00000000000000*Re(a0,0)

        """
        if ring is None:
            ring = self.parent()

        return ring({key: image for key, value in self._coefficients.items() if (image := f(value))})

    def __call__(self, values):
        r"""
        Return the value of this symbolic expression at ``values``.

        EXAMPLES::

            sage: from flatsurf.geometry.power_series import PowerSeriesCoefficientExpressionRing
            sage: R = PowerSeriesCoefficientExpressionRing(CC, gens=("Re(a0,?)", "Im(a0,?)"))
            sage: re = R.gen(("Re(a0,?)", 0))
            sage: im = R.gen(("Im(a0,?)", 0))

            sage: re({re: 1})
            1.00000000000000
            sage: re({re: 2, im: 1})
            2.00000000000000
            sage: (2 * re * im)({re: 3, im: 5})
            30.0000000000000

        """
        def evaluate(monomial):
            product = self.parent().base_ring().one()

            for variable in monomial:
                product *= values[self.parent().gen(variable)]

            return product

        return sum([
            coefficient * evaluate(monomial) for (monomial, coefficient) in self._coefficients.items()
            ])


class PowerSeriesCoefficientExpressionRing(UniqueRepresentation, CommutativeRing):
    r"""
    The ring of expressions in the symbolic coefficients of a multivariate
    power series.

    EXAMPLES::

        sage: from flatsurf.geometry.power_series import PowerSeriesCoefficientExpressionRing
        sage: R = PowerSeriesCoefficientExpressionRing(CC, gens=("Re(a0,?)", "Im(a0,?)"))
        sage: R
        Ring of Power Series Coefficients in Re(a0,0),…,Im(a0,0),… over Complex Field with 53 bits of precision

    TESTS::

        sage: R.has_coerce_map_from(CC)
        True
        sage: TestSuite(R).run()

    """

    def __init__(self, base_ring, gens, category=None):
        self._gens = gens

        from sage.categories.all import CommutativeRings

        CommutativeRing.__init__(self, base_ring, category=category or CommutativeRings(), normalize=False)
        self.register_coercion(base_ring)

    Element = PowerSeriesCoefficientExpression

    def _repr_(self):
        return f"Ring of Power Series Coefficients in {','.join(gen.replace('?', '0') + ',…' for gen in self._gens)} over {self.base_ring()}"

    def change_ring(self, ring):
        r"""
        Return this ring with the ring of constants replaced by ``ring``.

        EXAMPLES::

            sage: from flatsurf.geometry.power_series import PowerSeriesCoefficientExpressionRing
            sage: R = PowerSeriesCoefficientExpressionRing(CC, gens=("Re(a0,?)", "Im(a0,?)"))
            sage: R.change_ring(RR)
            Ring of Power Series Coefficients in Re(a0,0),…,Im(a0,0),… over Real Field with 53 bits of precision

        """
        return PowerSeriesCoefficientExpressionRing(ring, self._gens, category=self.category())

    def sum(self, summands):
        r"""
        Return the sum of ``summands``.

        This is an optimized version of the builtin `sum` that creates fewer
        temporary objects.

        EXAMPLES::

            sage: from flatsurf.geometry.power_series import PowerSeriesCoefficientExpressionRing
            sage: R = PowerSeriesCoefficientExpressionRing(CC, gens=("Re(a0,?)", "Im(a0,?)"))

            sage: R.sum([R.gen(0), R.gen(1)])
            Re(a0,0) + Im(a0,0)

        """
        # TODO: Add a benchmark to show that this is actually way faster when
        # there are lots of generators.

        summands = list(summands)

        if len(summands) == 0:
            return self.zero()

        if len(summands) == 1:
            return summands.pop()

        coefficients = dict(summands.pop()._coefficients)

        while summands:
            summand = summands.pop()
            for monomial, coefficient in summand._coefficients.items():
                assert coefficient
                if monomial not in coefficients:
                    coefficients[monomial] = coefficient
                else:
                    coefficients[monomial] += coefficient

                if not coefficients[monomial]:
                    del coefficients[monomial]

        return self(coefficients)

    def real_field(self):
        r"""
        Return a real base field corresponding to the complex base field of
        this ring.

        EXAMPLES::

            sage: from flatsurf.geometry.power_series import PowerSeriesCoefficientExpressionRing
            sage: R = PowerSeriesCoefficientExpressionRing(CC, gens=("Re(a0,?)", "Im(a0,?)"))

            sage: R.real_field()
            Real Field with 53 bits of precision

        When the base ring is not complex, this method is not functional::

            sage: R = PowerSeriesCoefficientExpressionRing(QQ, gens=("Re(a0,?)", "Im(a0,?)"))

            sage: R.real_field()
            Traceback (most recent call last):
            ...
            AttributeError: 'RationalField_with_category' object has no attribute 'prec'

        """
        from sage.all import RealField
        return RealField(self.base_ring().prec())

    def is_exact(self):
        r"""
        Return whether this ring is implementing exact arithmetic.

        EXAMPLES::

            sage: from flatsurf.geometry.power_series import PowerSeriesCoefficientExpressionRing
            sage: R = PowerSeriesCoefficientExpressionRing(QQ, gens=("a", "b"))
            sage: R.is_exact()
            True

        """
        return self.base_ring().is_exact()

    def _coerce_map_from_(self, other):
        if isinstance(other, PowerSeriesCoefficientExpressionRing):
            return self.base_ring().has_coerce_map_from(other.base_ring())

    def _element_constructor_(self, x):
        r"""
        Return an element of this ring built from ``x``.

        EXAMPLES::

            sage: from flatsurf.geometry.power_series import PowerSeriesCoefficientExpressionRing
            sage: R = PowerSeriesCoefficientExpressionRing(QQ, gens=("a?", "b?"))
            sage: R(1)
            1
            sage: R({(1,): 3})
            3*b0

        """
        if isinstance(x, PowerSeriesCoefficientExpression):
            return x.map_coefficients(self.base_ring(), ring=self)

        if isinstance(x, dict):
            return self.element_class(self, {
                tuple(sorted(monomial)): self.base_ring()(coefficient) for (monomial, coefficient) in x.items() if coefficient
            })

        from sage.all import parent
        if parent(x) is self.base_ring():
            if not x:
                return self.element_class(self, {})
            return self.element_class(self, {(): x})

        raise TypeError(f"cannot create a symbolic expression from this {type(x)}")

    @cached_method
    def gen(self, gen):
        r"""
        Return the generator identified by ``gen``.

        INPUT:

        - ``gen`` -- a tuple (name, degree) or an integer.

        EXAMPLES::

            sage: from flatsurf.geometry.power_series import PowerSeriesCoefficientExpressionRing
            sage: R = PowerSeriesCoefficientExpressionRing(QQ, gens=("a?", "b?"))
            sage: R.gen(("a?", 0))
            a0

        SageMath wants us to also order the generators and return them when
        calling ``gen`` with an integer argument::

            sage: R.gen(0)
            a0
            sage: R.gen(1)
            b0
            sage: R.gen(2)
            a1
            sage: R.gen(3)
            b1

        """
        if isinstance(gen, tuple):
            if len(gen) == 2:
                name, degree = gen
                if name not in self._gens:
                    raise ValueError(f"name must be one of {self._gens} not {name}")
                gen = degree * len(self._gens) + self._gens.index(name)

        from sage.all import parent, ZZ
        if parent(gen) == ZZ:
            gen = int(gen)

        if isinstance(gen, int):
            return self.element_class(self, {(gen,): self.base().one()})

        raise TypeError("gen must be an integer or a tuple (name, degree)")

    def ngens(self):
        r"""
        Return the number of generators of this ring.

        Since there are infinitely many generators, this method is not
        implemented (SageMath does not accept us returning +infinity here.)

        EXAMPLES::

            sage: from flatsurf.geometry.power_series import PowerSeriesCoefficientExpressionRing
            sage: R = PowerSeriesCoefficientExpressionRing(QQ, gens=("a?", "b?"))
            sage: R.ngens()
            Traceback (most recent call last):
            ...
            NotImplementedError

        """
        raise NotImplementedError
