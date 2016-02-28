#include <math.h>
#include <SU3NCSMUtils/factorial.h>
#include <SU3NCSMUtils/clebschGordan.h>

#ifdef CWig9jLookUpTable_PROFILE
int CWig9jLookUpTable::nbfetch=0;
int CWig9jLookUpTable::nbadd=0;
#endif

template< class T > inline const T& max(const T& x, const T& y) { return ( x > y ? x : y ); }
template< class T > inline const T& min(const T& x, const T& y) { return ( x < y ? x : y ); }

/*
 * Clebsch Gordan coefficients
 * -------------------------------
 * Update: C. Bahri (U of T, 2/99)
 *
 * Ref: B.L. Van der Waerden, Group Theoretical Methods in Quantum Mechanics
        G. Racah, Phys. Rev. 62 (1942) 438
        D.A. Varshalovich et al, Quantum Theory of Angular Momentum
            (World Scientific 1988) p. 238
 */

double clebschGordan(int a2, int al2, int b2, int bt2,
                     int c2, int gm2)
{
    double result = 0., numerator, denominator = 0.;

    if (isAllowed(a2, b2, c2) && gm2 == al2 + bt2
        && !isOdd(a2+al2)     && !isOdd(b2+bt2))
    {
        int zmin, zmax;

        zmin = max(max(0, -(c2 - b2 + al2)), -(c2 - a2 - bt2));
        zmax = min(min(a2 + b2 - c2, a2 - al2), b2 + bt2);

        numerator = 0.5 * (  logDelta2(a2, b2, c2) + log(double (c2 + 1))
                           + logFact((a2 + al2)/2) + logFact((a2 - al2)/2)
                           + logFact((b2 + bt2)/2) + logFact((b2 - bt2)/2)
                           + logFact((c2 + gm2)/2) + logFact((c2 - gm2)/2));

        for (int z2 = zmin; z2 <= zmax; z2 += 2)
        {
            denominator =   logFact(z2/2)
                          + logFact((a2 + b2 - c2 - z2)/2)
                          + logFact((a2 - al2 - z2)/2)
                          + logFact((b2 + bt2 - z2)/2)
                          + logFact((c2 - b2 + al2 + z2)/2)
                          + logFact((c2 - a2 - bt2 + z2)/2);

            if (isOdd(z2/2))
                result -= exp(numerator - denominator);
            else
                result += exp(numerator - denominator);
        }
    }
    return (result);
} /* clebschGordan */

/*
 * Racah coefficients
 * -------------------------------
 * Update: C. Bahri (U of T, 2/99)
 *
 * Ref: G. Racah, Phys. Rev. 62 (1942) 438
        D.A. Varshalovich et al, Quantum Theory of Angular Momentum
            (World Scientific 1988) p. 293
 */

double racah(int a2, int b2, int e2, int d2, int c2, int f2)
{
    double result = 0., numerator, denominator = 0.;

    if (   isAllowed(a2, b2, c2) && isAllowed(c2, d2, e2)
        && isAllowed(a2, e2, f2) && isAllowed(b2, d2, f2))
    {
        int zmin, zmax;
//      register int

        zmin = max(max(0, -(-a2 + c2 - d2 + f2)), -(-b2 + c2 - e2 + f2));
        zmax = min(min(min(min(a2 + b2 + d2 + e2 + 2,
               a2 + b2 - c2), -c2 + d2 + e2), a2 + e2 - f2), b2 + d2 - f2);

        numerator = 0.5 * (  logDelta2(a2, b2, c2) + logDelta2(c2, d2, e2)
                           + logDelta2(a2, e2, f2) + logDelta2(b2, d2, f2));

        for (int z2 = zmin; z2 <= zmax; z2 += 2)
        {
            denominator = logFact((a2 + b2 + d2 + e2 + 2 - z2)/2)
                          - (  logFact(z2/2)
                             + logFact(( a2 + b2 - c2 - z2)/2)
                             + logFact((-c2 + d2 + e2 - z2)/2)
                             + logFact(( a2 + e2 - f2 - z2)/2)
                             + logFact(( b2 + d2 - f2 - z2)/2)
                             + logFact((-a2 + c2 - d2 + f2 + z2)/2)
                             + logFact((-b2 + c2 - e2 + f2 + z2)/2));

            if (isOdd(z2/2))
                result -= exp(numerator + denominator);
            else
                result += exp(numerator + denominator);
        }
    }
    return (result);
} /* racah */

/*
 * Wigner symbols
 * -------------------------------
 * Update: C. Bahri (U of T, 2/99)
 *
 * Ref:
 */

double wigner3jm(int a2,  int b2,  int c2,
                 int al2, int bt2, int gm2)
{
    double result = 0., numerator, denominator = 0.;

    if (isAllowed(a2, b2, c2) && al2 + bt2 + gm2 == 0
        && !isOdd(a2+al2)     && !isOdd(b2+bt2))
    {
        int zmin, zmax;

        zmin = max(max(0, -(c2 - b2 + al2)), -(c2 - a2 - bt2));
        zmax = min(min(a2 + b2 - c2, a2 - al2), b2 + bt2);

        numerator = 0.5 * (  logDelta2(a2, b2, c2)
                           + logFact((a2 + al2)/2) + logFact((a2 - al2)/2)
                           + logFact((b2 + bt2)/2) + logFact((b2 - bt2)/2)
                           + logFact((c2 + gm2)/2) + logFact((c2 - gm2)/2));

        for (int z2 = zmin; z2 <= zmax; z2 += 2)
        {
            denominator =   logFact(z2/2)
                          + logFact((a2 + b2 - c2 - z2)/2)
                          + logFact((a2 - al2 - z2)/2)
                          + logFact((b2 + bt2 - z2)/2)
                          + logFact((c2 - b2 + al2 + z2)/2)
                          + logFact((c2 - a2 - bt2 + z2)/2);

            if (isOdd2(z2))
                result -= exp(numerator - denominator);
            else
                result += exp(numerator - denominator);
        }
        if (isOdd2(a2 - b2 + gm2)) result = -result;
    }
    return (result);
} /* wigner3jm */


double wigner6j(int a2, int b2, int c2,
                int d2, int e2, int f2)
{
    double result = 0., numerator, denominator = 0.;

    if (   isAllowed(a2, b2, c2) && isAllowed(c2, d2, e2)
        && isAllowed(a2, e2, f2) && isAllowed(b2, d2, f2))
    {
        int zmin, zmax;
        register int abc, cde, aef, bdf, abde, acdf, bcef;

        zmin = max(max(max(abc=(a2 + b2 + c2)/2, cde=(c2 + d2 + e2)/2),
               aef=(a2 + e2 + f2)/2), bdf=(b2 + d2 + f2)/2);
        zmax = min(min(abde=(a2 + b2 + d2 + e2)/2, acdf=(a2 + c2 + d2 + f2)/2),
               bcef=(b2 + c2 + e2 + f2)/2);

        numerator = 0.5 * (  logDelta2(a2, b2, c2) + logDelta2(c2, d2, e2)
                           + logDelta2(a2, e2, f2) + logDelta2(b2, d2, f2));

        for (int z = zmin; z <= zmax; ++z)
        {
            denominator = logFact(z + 1)
                          - (  logFact(z - abc)
                             + logFact(z - cde)
                             + logFact(z - aef)
                             + logFact(z - bdf)
                             + logFact(abde - z)
                             + logFact(acdf - z)
                             + logFact(bcef - z));

            if (isOdd(z))
                result -= exp(numerator + denominator);
            else
                result += exp(numerator + denominator);
        }
    }
    return (result);
} /* wigner6j */


double wigner9j(int a2, int b2, int c2,
                int d2, int e2, int f2,
                int g2, int h2, int j2)
{
    double result = 0.;

    if (   isAllowed(a2, b2, c2) && isAllowed(d2, e2, f2) 
        && isAllowed(g2, h2, j2) && isAllowed(a2, d2, g2)
        && isAllowed(b2, e2, h2) && isAllowed(c2, f2, j2))
    {
        int zmin, zmax;

        zmin = max(max(abs(a2 - j2), abs(b2 - f2)), abs(d2 - h2));
        zmax = min(min(    a2 + j2 ,     b2 + f2 ),     d2 + h2 );

        for (int z2 = zmin; z2 <= zmax; z2 += 2)
        {
            result += double (z2 + 1) * racah(a2, j2, d2, h2, z2, g2)
                                      * racah(a2, j2, b2, f2, z2, c2)
                                      * racah(b2, f2, h2, d2, z2, e2);
        }
    }
    return (result);
} /* wigner9j */

/*
 * Unitary coefficients
 * -------------------------------
 * Update: C. Bahri (U of T, 3/99)
 *
 * Ref:
 */

double racahU(int a2, int b2, int e2, int d2, int c2, int f2)
{
    double result = 0., numerator, denominator = 0.;

    if (isAllowed(a2, b2, c2) && isAllowed(c2, d2, e2)
        && isAllowed(a2, e2, f2) && isAllowed(b2, d2, f2))
    {
        int zmin, zmax;
//      register int

        zmin = max(max(0, -(-a2 + c2 - d2 + f2)), -(-b2 + c2 - e2 + f2));
        zmax = min(min(min(min(a2 + b2 + d2 + e2 + 2,
               a2 + b2 - c2), -c2 + d2 + e2), a2 + e2 - f2), b2 + d2 - f2);

        numerator = 0.5 * (  logDelta2(a2, b2, c2) + logDelta2(c2, d2, e2)
                           + logDelta2(a2, e2, f2) + logDelta2(b2, d2, f2)
                           + log (double (c2 + 1)) + log (double (f2 + 1)));

        for (int z2 = zmin; z2 <= zmax; z2 += 2)
        {
            denominator = logFact((a2 + b2 + d2 + e2 + 2 - z2)/2)
                          - (  logFact(z2/2)
                             + logFact(( a2 + b2 - c2 - z2)/2)
                             + logFact((-c2 + d2 + e2 - z2)/2)
                             + logFact(( a2 + e2 - f2 - z2)/2)
                             + logFact(( b2 + d2 - f2 - z2)/2)
                             + logFact((-a2 + c2 - d2 + f2 + z2)/2)
                             + logFact((-b2 + c2 - e2 + f2 + z2)/2));

            if (isOdd2(z2))
                result -= exp(numerator + denominator);
            else
                result += exp(numerator + denominator);
        }
    }
    return (result);
} /* racahU */


double unitary9j(int a2, int b2, int c2,
                 int d2, int e2, int f2,
                 int g2, int h2, int j2)
{
    double result = 0.;

    if (isAllowed(a2, b2, c2) && isAllowed(d2, e2, f2) 
        && isAllowed(g2, h2, j2) && isAllowed(a2, d2, g2)
        && isAllowed(b2, e2, h2) && isAllowed(c2, f2, j2))
    {
        int zmin, zmax;

        zmin = max(max(abs(a2 - j2), abs(b2 - f2)), abs(d2 - h2));
        zmax = min(min(    a2 + j2 ,     b2 + f2 ),     d2 + h2 );

        for (int z2 = zmin; z2 <= zmax; z2 += 2)
        {
            result += double (z2 + 1) * racah(a2, j2, d2, h2, z2, g2)
                                      * racah(a2, j2, b2, f2, z2, c2)
                                      * racah(b2, f2, h2, d2, z2, e2);
        }
        result *= sqrt(double ((c2 + 1)*(f2 + 1)*(g2 + 1)*(h2 + 1)));
    }
    return (result);
} /* unitary9j */

