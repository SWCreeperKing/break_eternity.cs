using System.Globalization;

namespace Break_Eternity;

public class Eternity
{
    public enum Format
    {
        Mixed = 0,
        _ = Mixed
    }

    public static Format DecFormat { get; set; } = Format.Mixed;

    private static readonly double MaxSignificantDigits = 17;
    private static readonly double ExponentLimit = 9e15;
    private static readonly double LayerDown = 15.954242509439325;
    private static readonly double FirstNegativeLayer = 1.1111111111111112e-16;
    private static readonly double NumberExponentMin = -324;
    private static readonly double NumberExponentMax = 308;
    private static readonly double MaxEsInARow = 5;
    private static readonly double Omega = 0.56714329040978387299997;
    private static readonly bool IgnoreCommas = true;
    private static readonly bool CommasAreDecPoints = false;

    public static readonly Eternity Zero = FC_NN(0, 0, 0);
    public static readonly Eternity One = FC_NN(1, 0, 1);
    public static readonly Eternity NegativeOne = FC_NN(-1, 0, 1);
    public static readonly Eternity Two = FC_NN(1, 0, 2);
    public static readonly Eternity Ten = FC_NN(1, 0, 10);
    public static readonly Eternity NaN = FC_NN(int.MaxValue, long.MaxValue, double.NaN);
    public static readonly Eternity Infinity = FC_NN(1, long.MaxValue, double.PositiveInfinity);
    public static readonly Eternity NegativeInfinity = FC_NN(-1, long.MinValue, double.NegativeInfinity);
    public static readonly Eternity NumberMax = FC(1, 0, double.MaxValue);
    public static readonly Eternity NumberMin = FC(1, 0, double.MinValue);
    public static readonly Eternity MathE = new(Math.E);

    private static readonly double[] CriticalHeaders = { 2, 2.718281828459045, 3, 4, 5, 6, 7, 8, 9, 10 };

    private static readonly double[,] CriticalTetrValues =
    {
        {
            1,
            1.0891168053867777,
            1.1789745164521264,
            1.2701428397304229,
            1.3632066654400328,
            1.4587804913784246,
            1.557523817412741,
            1.660158301473385,
            1.767487542936873,
            1.8804205225512542,
            2
        },
        {
            1,
            1.1121114330934,
            1.23103892493161,
            1.35838369631113,
            1.49605193039935,
            1.64635423375119,
            1.81213853570186,
            1.99697132461829,
            2.20538955455724,
            2.44325744833852,
            2.718281828459045
        },
        {
            1,
            1.1187738849693603,
            1.2464963939368214,
            1.38527004705667,
            1.5376664685821402,
            1.7068895236551784,
            1.897001227148399,
            2.1132403089001035,
            2.362480153784171,
            2.6539010333870774,
            3
        },
        {
            1,
            1.1367350847096405,
            1.2889510672956703,
            1.4606478703324786,
            1.6570295196661111,
            1.8850062585672889,
            2.1539465047453485,
            2.476829779693097,
            2.872061932789197,
            3.3664204535587183,
            4
        },
        {
            1,
            1.1494592900767588,
            1.319708228183931,
            1.5166291280087583,
            1.748171114438024,
            2.0253263297298045,
            2.3636668498288547,
            2.7858359149579424,
            3.3257226212448145,
            4.035730287722532,
            5
        },
        {
            1,
            1.159225940787673,
            1.343712473580932,
            1.5611293155111927,
            1.8221199554561318,
            2.14183924486326,
            2.542468319282638,
            3.0574682501653316,
            3.7390572020926873,
            4.6719550537360774,
            6
        },
        {
            1,
            1.1670905356972596,
            1.3632807444991446,
            1.5979222279405536,
            1.8842640123816674,
            2.2416069644878687,
            2.69893426559423,
            3.3012632110403577,
            4.121250340630164,
            5.281493033448316,
            7
        },
        {
            1,
            1.1736630594087796,
            1.379783782386201,
            1.6292821855668218,
            1.9378971836180754,
            2.3289975651071977,
            2.8384347394720835,
            3.5232708454565906,
            4.478242031114584,
            5.868592169644505,
            8
        },
        {
            1,
            1.1793017514670474,
            1.394054150657457,
            1.65664127441059,
            1.985170999970283,
            2.4069682290577457,
            2.9647310119960752,
            3.7278665320924946,
            4.814462547283592,
            6.436522247411611,
            9
        },
        {
            1,
            1.18422737399915,
            1.4066113788546144,
            1.680911177655277,
            2.027492094355525,
            2.4775152854601967,
            3.080455730250329,
            3.918234505962507,
            5.1332705696484595,
            6.9878696918072905,
            10
        }
    };

    private static readonly double[,] CriticalSLogValues =
    {
        {
            -1,
            -0.9194161097107025,
            -0.8335625019330468,
            -0.7425599821143978,
            -0.6466611521029437,
            -0.5462617907227869,
            -0.4419033816638769,
            -0.3342645487554494,
            -0.224140440909962,
            -0.11241087890006762,
            0
        },
        {
            -1,
            -0.90603157029014,
            -0.80786507256596,
            -0.7064666939634,
            -0.60294836853664,
            -0.49849837513117,
            -0.39430303318768,
            -0.29147201034755,
            -0.19097820800866,
            -0.09361896280296,
            0
        },
        {
            -1,
            -0.9021579584316141,
            -0.8005762598234203,
            -0.6964780623319391,
            -0.5911906810998454,
            -0.486050182576545,
            -0.3823089430815083,
            -0.28106046722897615,
            -0.1831906535795894,
            -0.08935809204418144,
            0
        },
        {
            -1,
            -0.8917227442365535,
            -0.781258746326964,
            -0.6705130326902455,
            -0.5612813129406509,
            -0.4551067709033134,
            -0.35319256652135966,
            -0.2563741554088552,
            -0.1651412821106526,
            -0.0796919581982668,
            0
        },
        {
            -1,
            -0.8843387974366064,
            -0.7678744063886243,
            -0.6529563724510552,
            -0.5415870994657841,
            -0.4352842206588936,
            -0.33504449124791424,
            -0.24138853420685147,
            -0.15445285440944467,
            -0.07409659641336663,
            0
        },
        {
            -1,
            -0.8786709358426346,
            -0.7577735191184886,
            -0.6399546189952064,
            -0.527284921869926,
            -0.4211627631006314,
            -0.3223479611761232,
            -0.23107655627789858,
            -0.1472057700818259,
            -0.07035171210706326,
            0
        },
        {
            -1,
            -0.8740862815291583,
            -0.7497032990976209,
            -0.6297119746181752,
            -0.5161838335958787,
            -0.41036238255751956,
            -0.31277212146489963,
            -0.2233976621705518,
            -0.1418697367979619,
            -0.06762117662323441,
            0
        },
        {
            -1,
            -0.8702632331800649,
            -0.7430366914122081,
            -0.6213373075161548,
            -0.5072025698095242,
            -0.40171437727184167,
            -0.30517930701410456,
            -0.21736343968190863,
            -0.137710238299109,
            -0.06550774483471955,
            0
        },
        {
            -1,
            -0.8670016295947213,
            -0.7373984232432306,
            -0.6143173985094293,
            -0.49973884395492807,
            -0.394584953527678,
            -0.2989649949848695,
            -0.21245647317021688,
            -0.13434688362382652,
            -0.0638072667348083,
            0
        },
        {
            -1,
            -0.8641642839543857,
            -0.732534623168535,
            -0.6083127477059322,
            -0.4934049257184696,
            -0.3885773075899922,
            -0.29376029055315767,
            -0.2083678561173622,
            -0.13155653399373268,
            -0.062401588652553186,
            0
        }
    };

    private static readonly List<double> PowersOf10 = PowerOf10();
    private static double GetPowerOf10(int index) => PowersOf10[index + 323];

    private static List<double> PowerOf10()
    {
        List<double> powers = new();
        for (var i = (int) (NumberExponentMin + 1); i <= NumberExponentMax; i++)
        {
            powers.Add(double.Parse($"1e{i}", System.Globalization.NumberStyles.Float));
        }

        return powers;
    }

    public int S
    {
        get => _sign;
        set
        {
            if (value == 0)
            {
                _sign = 0;
                _layer = 0;
                _mag = 0;
            }
            else _sign = value;
        }
    }

    public double M
    {
        get
        {
            if (_sign == 0) return 0;
            switch (_layer)
            {
                case 0:
                    var exp = Math.Floor(Math.Log10(_mag));
                    double man;
                    if (_mag == 5e-324) man = 5;
                    else man = _mag / GetPowerOf10((int) exp);
                    return _sign * man;
                case 1:
                    var residue = _mag - Math.Floor(_mag);
                    return _sign * Math.Pow(10, residue);
                default:
                    return _sign;
            }
        }
        set
        {
            if (_layer <= 2) FromMantissaExponent(value, E);
            else
            {
                _sign = Math.Sign(value);
                if (_sign != 0) return;
                _layer = 0;
                E = 0;
            }
        }
    }

    public double E
    {
        get
        {
            if (_sign == 0) return 0;
            return _layer switch
            {
                0 => Math.Floor(Math.Log10(_mag)),
                1 => Math.Floor(_mag),
                2 => Math.Floor(Math.Sign(_mag) * Math.Pow(10, Math.Abs(_mag))),
                _ => _mag * double.PositiveInfinity
            };
        }
        set => FromMantissaExponent(M, value);
    }

    private int _sign;
    private double _mag;
    private long _layer;

    public Eternity()
    {
    }

    public Eternity(int sign, long layer, double magnitude, bool normalize = true)
    {
        if (normalize) FromComponents(sign, layer, magnitude);
        else FromComponentsNoNormalize(sign, layer, magnitude);
    }

    public Eternity(double mantissa, double exponent) => FromMantissaExponent(mantissa, exponent);
    public Eternity(string s) => FromString(s);
    public Eternity(double d) => DE(d);
    public Eternity(decimal d) => DE((double) d);
    public Eternity(long d) => DE(d);
    public Eternity(ulong d) => DE(d);

    public static implicit operator Eternity(string s) => new(s);
    public static implicit operator Eternity(decimal d) => new(d);
    public static implicit operator Eternity(double d) => new(d);
    public static implicit operator Eternity(float d) => new(d);
    public static implicit operator Eternity(sbyte d) => new(d);
    public static implicit operator Eternity(byte d) => new(d);
    public static implicit operator Eternity(short d) => new(d);
    public static implicit operator Eternity(ushort d) => new(d);
    public static implicit operator Eternity(int d) => new(d);
    public static implicit operator Eternity(uint d) => new(d);
    public static implicit operator Eternity(long d) => new(d);
    public static implicit operator Eternity(ulong d) => new(d);

    public static explicit operator decimal(Eternity d) => (decimal) d.ToNumber();
    public static explicit operator double(Eternity d) => d.ToNumber();
    public static explicit operator float(Eternity d) => (float) d.ToNumber();
    public static explicit operator sbyte(Eternity d) => (sbyte) d.ToNumber();
    public static explicit operator byte(Eternity d) => (byte) d.ToNumber();
    public static explicit operator short(Eternity d) => (short) d.ToNumber();
    public static explicit operator ushort(Eternity d) => (ushort) d.ToNumber();
    public static explicit operator int(Eternity d) => (int) d.ToNumber();
    public static explicit operator uint(Eternity d) => (uint) d.ToNumber();
    public static explicit operator long(Eternity d) => (long) d.ToNumber();
    public static explicit operator ulong(Eternity d) => (ulong) d.ToNumber();

    public static bool operator ==(Eternity d1, Eternity d2) => Eq(d1, d2);
    public static bool operator !=(Eternity d1, Eternity d2) => Neq(d1, d2);
    public static bool operator >(Eternity d1, Eternity d2) => Gt(d1, d2);
    public static bool operator <(Eternity d1, Eternity d2) => Lt(d1, d2);
    public static bool operator >=(Eternity d1, Eternity d2) => Gte(d1, d2);
    public static bool operator <=(Eternity d1, Eternity d2) => Lte(d1, d2);

    public static Eternity operator -(Eternity d) => Neg(d);
    public static Eternity operator +(Eternity d1, Eternity d2) => d1.Add(d2);
    public static Eternity operator -(Eternity d1, Eternity d2) => d1.Sub(d2);
    public static Eternity operator *(Eternity d1, Eternity d2) => d1.Mul(d2);
    public static Eternity operator /(Eternity d1, Eternity d2) => d1.Div(d2);
    public static Eternity operator ++(Eternity d) => Add(d, 1);
    public static Eternity operator --(Eternity d) => Sub(d, 1);

    public static Eternity D(Eternity d) => new Eternity().FromDec(d);
    public static Eternity N(double d) => new Eternity().FromNumber(d);
    public static Eternity FC(int sign, long layer, double mag) => new Eternity().FromComponents(sign, layer, mag);

    public static Eternity FC_NN(int sign, long layer, double mag) =>
        new Eternity().FromComponentsNoNormalize(sign, layer, mag);

    public static Eternity ME(double man, double exp) => new Eternity().FromMantissaExponent(man, exp);
    public static Eternity Abs(Eternity d) => D(d).Abs();
    public static Eternity Neg(Eternity d) => D(d).Neg();
    public static Eternity Round(Eternity d) => D(d).Round();
    public static Eternity Floor(Eternity d) => D(d).Floor();
    public static Eternity Ceil(Eternity d) => D(d).Ceil();
    public static Eternity Trunc(Eternity d) => D(d).Trunc();
    public static Eternity Add(Eternity d1, Eternity d2) => d1.Add(d2);
    public static Eternity Sub(Eternity d1, Eternity d2) => D(d1).Sub(d2);
    public static Eternity Mul(Eternity d1, Eternity d2) => D(d1).Mul(d2);
    public static Eternity Div(Eternity d1, Eternity d2) => D(d1).Div(d2);
    public static Eternity Max(Eternity d1, Eternity d2) => D(d1).Max(d2);
    public static Eternity Min(Eternity d1, Eternity d2) => D(d1).Min(d2);
    public static Eternity AbsMin(Eternity d1, Eternity d2) => D(d1).AbsMin(d2);
    public static Eternity AbsMax(Eternity d1, Eternity d2) => D(d1).AbsMax(d2);
    public static Eternity Clamp(Eternity d1, Eternity d2, Eternity d3) => D(d1).Clamp(d2, d3);
    public static Eternity AbsLog10(Eternity d) => D(d).AbsLog10();
    public static Eternity Log10(Eternity d) => D(d).Log10();
    public static Eternity Log(Eternity d1, Eternity d2) => D(d1).Log(d2);
    public static Eternity Log2(Eternity d) => D(d).Log2();
    public static Eternity Ln(Eternity d) => D(d).Ln();
    public static Eternity Pow(Eternity d1, Eternity d2) => D(d1).Pow(d2);
    public static Eternity Pow10(Eternity d) => D(d).Pow10();
    public static Eternity Root(Eternity d1, Eternity d2) => D(d1).Root(d2);
    public static Eternity Factorial(Eternity d) => D(d).Factorial();
    public static Eternity Gamma(Eternity d) => D(d).Gamma();
    public static Eternity LnGamma(Eternity d) => D(d).LnGamma();
    public static Eternity Exp(Eternity d) => D(d).Exp();
    public static Eternity Sq(Eternity d) => D(d).Sq();
    public static Eternity Sqrt(Eternity d) => D(d).Sqrt();
    public static Eternity Cube(Eternity d) => D(d).Cube();
    public static Eternity Cbrt(Eternity d) => D(d).Cbrt();

    public static Eternity Tetrate(Eternity d, double h = 2, Eternity? p = null) =>
        D(d).Tetrate(h, p ?? FC_NN(1, 0, 1));

    public static Eternity IteratedLog(Eternity d, Eternity? b = null, double t = 1) =>
        D(d).IteratedLog(b ?? FC_NN(1, 0, 10), t);

    public static Eternity LayerAdd10(Eternity d1, double diff) => D(d1).LayerAdd10(diff);

    public static Eternity LayerAdd(Eternity d, double diff, Eternity? b = null) =>
        D(d).LayerAdd(diff, b ?? FC_NN(1, 0, 10));

    public static Eternity SLog(Eternity d, Eternity? b = null) => SLog(d, b ?? FC_NN(1, 0, 10));
    public static Eternity LambertW(Eternity d) => D(d).LambertW();

    public static Eternity DLambertW(Eternity z, double tol = 1e-10)
    {
        Eternity w, ew, wewz, wn;
        if (!double.IsFinite(z._mag)) return z;
        if (z == Zero) return z;
        if (z == One) return D(Omega);
        var absz = Abs(z);
        w = Ln(z);
        for (var i = 0; i < 100; i++)
        {
            ew = Exp(-w);
            wewz = w.Sub(z.Mul(ew));
            wn = w.Sub(wewz.Div(w.Add(1).Sub(w.Add(2).Mul(wewz).Div(Mul(2, w).Add(2)))));
            if (Abs(wn.Sub(w)).Lt(Abs(wn).Mul(tol))) return wn;
            w = wn;
        }

        return NaN;
    }

    public static Eternity SSqrt(Eternity d) => D(d).SSqrt();

    public static Eternity Pentate(Eternity d, double h = 2, Eternity? p = null) =>
        D(d).Pentate(h, p ?? FC_NN(1, 0, 1));

    public static int Sign(Eternity d) => d._sign;
    public static int Cmp(Eternity d1, Eternity d2) => D(d1).Cmp(d2);
    public static int CmpAbs(Eternity d1, Eternity d2) => D(d1).CmpAbs(d2);
    public static bool IsNaN(Eternity d) => double.IsNaN(d._mag);
    public static bool IsFinite(Eternity d) => double.IsFinite(d._mag);
    public static bool Eq(Eternity d1, Eternity d2) => D(d1).Eq(d2);
    public static bool Neq(Eternity d1, Eternity d2) => D(d1).Neq(d2);
    public static bool Lt(Eternity d1, Eternity d2) => D(d1).Lt(d2);
    public static bool Lte(Eternity d1, Eternity d2) => D(d1).Lte(d2);
    public static bool Gt(Eternity d1, Eternity d2) => D(d1).Gt(d2);
    public static bool Gte(Eternity d1, Eternity d2) => D(d1).Gte(d2);
    public static string FormatDec(Eternity d, int? p = null, Format? f = null) => d.FormatDec(p, f);

    public static double SLogCritical(double b, double height)
    {
        return b > 10 ? height - 1 : CriticalSection(b, height, CriticalSLogValues);
    }

    public static double TetrateCritical(double b, double height) => CriticalSection(b, height, CriticalTetrValues);

    public static double CriticalSection(double b, double height, double[,] grid)
    {
        height *= 10;
        if (height < 0) height = 0;
        if (height > 10) height = 10;
        if (b < 2) b = 2;
        if (b > 10) b = 10;

        double lower = 0;
        double upper = 0;
        for (var i = 0; i < CriticalHeaders.Length; i++)
        {
            if (CriticalHeaders[i] == b)
            {
                lower = grid[i, (int) Math.Floor(height)];
                upper = grid[i, (int) Math.Ceiling(height)];
                break;
            }

            if (!(CriticalHeaders[i] < b) || !(CriticalHeaders[i + 1] > b)) continue;
            var basefrac = (b - CriticalHeaders[i]) / (CriticalHeaders[i + 1] - CriticalHeaders[i]);
            lower = grid[i, (int) Math.Floor(height)] * (1 - basefrac) +
                    grid[i + 1, (int) Math.Floor(height)] * basefrac;
            upper = grid[i, (int) Math.Ceiling(height)] * (1 - basefrac) +
                    grid[i + 1, (int) Math.Ceiling(height)] * basefrac;
            break;
        }

        var frac = height - Math.Floor(height);
        var result = lower * (1 - frac) + upper * frac;
        return result;
    }

    public static double DecPlaces(double value, int places)
    {
        var len = places + 1;
        var numDigits = Math.Ceiling(Math.Log10(Math.Abs(value)));
        var rounded = Math.Round(value * Math.Pow(10, len - numDigits)) * Math.Pow(10, numDigits - len);
        return double.Parse(rounded.ToString($"N{(int) Math.Max(len - numDigits, 0)}"));
    }

    public static double FMagLog10(double n) => Math.Sign(n) * Math.Log10(Math.Abs(n));

    public static double FGamma(double n)
    {
        if (!double.IsFinite(n)) return n;
        if (n < -50)
        {
            if (n == Math.Truncate(n)) return double.NegativeInfinity;
            return 0;
        }

        double scal1 = 1;
        while (n < 10)
        {
            scal1 *= n;
            n++;
        }

        n -= 1;
        var l = 0.9189385332046727;
        l += (n + 0.5) * Math.Log(n);
        l -= n;
        var n2 = n * n;
        var np = n;
        l += 1 / (12 * np);
        np *= n2;
        l += 1 / (360 * np);
        np *= n2;
        l += 1 / (1260 * np);
        np *= n2;
        l += 1 / (1680 * np);
        np *= n2;
        l += 1 / (1188 * np);
        np *= n2;
        l += 691 / (360360 * np);
        np *= n2;
        l += 7 / (1092 * np);
        np *= n2;
        l += 3617 / (122400 * np);
        return Math.Exp(l) / scal1;
    }

    public static double FLambertW(double z, double tol = 1e-10)
    {
        double w;
        double wn;
        if (!double.IsFinite(z)) return z;

        switch (z)
        {
            case 0:
                return z;
            case 1:
                return Omega;
            case < 10:
                w = 0;
                break;
            default:
                w = Math.Log(z) - Math.Log(Math.Log(z));
                break;
        }

        for (var i = 0; i < 100; i++)
        {
            wn = (z * Math.Exp(-w) + w * w) / (w + 1);
            if (Math.Abs(wn - w) < tol * Math.Abs(wn)) return wn;
            w = wn;
        }

        return double.NaN;
    }

    public Eternity DE(double d)
    {
        var n = new Eternity().FromNumber(d);
        FromComponentsNoNormalize(n._sign, n._layer, n._mag);
        return this;
    }

    public Eternity Normalize()
    {
        if (_sign == 0 || (_mag == 0 && _layer == 0))
        {
            _sign = 0;
            _mag = 0;
            _layer = 0;
            return this;
        }

        if (_layer == 0 && _mag < 0)
        {
            _mag = -_mag;
            _sign = -_sign;
        }

        if (_layer == 0 && _mag < FirstNegativeLayer)
        {
            _layer += 1;
            _mag = Math.Log10(_mag);
            return this;
        }

        var absmag = Math.Abs(_mag);
        double signmag = Math.Sign(_mag);
        if (absmag >= ExponentLimit)
        {
            _layer += 1;
            _mag = signmag * Math.Log10(absmag);
            return this;
        }

        while (absmag < LayerDown && _layer > 0)
        {
            _layer -= 1;
            if (_layer == 0)
                _mag = Math.Pow(10, _mag);
            else
            {
                _mag = signmag * Math.Pow(10, absmag);
                absmag = Math.Abs(_mag);
                signmag = Math.Sign(_mag);
            }
        }

        if (_layer != 0) return this;

        switch (_mag)
        {
            case < 0:
                _mag = -_mag;
                _sign = -_sign;
                break;
            case 0:
                _sign = 0;
                break;
        }

        return this;
    }

    public Eternity FromComponents(int sign, long layer, double mag)
    {
        _sign = sign;
        _layer = layer;
        _mag = mag;

        Normalize();
        return this;
    }

    public Eternity FromComponentsNoNormalize(int sign, long layer, double mag)
    {
        _sign = sign;
        _layer = layer;
        _mag = mag;
        return this;
    }

    public Eternity FromMantissaExponent(double man, double exp)
    {
        _layer = 1;
        _sign = Math.Sign(man);
        _mag = exp + Math.Log10(Math.Abs(man));

        Normalize();
        return this;
    }

    public Eternity FromDec(Eternity d)
    {
        _sign = d._sign;
        _layer = d._layer;
        _mag = d._mag;
        return this;
    }

    public Eternity FromNumber(double d)
    {
        _mag = Math.Abs(d);
        _sign = Math.Sign(d);
        _layer = 0;
        Normalize();
        return this;
    }

    public Eternity FromString(string value)
    {
        if (IgnoreCommas) value = value.Replace(",", "");
        else if (CommasAreDecPoints) value = value.Replace(',', '.');

        var pentationparts = value.Split("^^^");
        if (pentationparts.Length == 2)
        {
            var b = double.Parse(pentationparts[0]);
            var h = double.Parse(pentationparts[1]);
            var heightparts = pentationparts[1].Split(";");
            double payload = 1;

            if (heightparts.Length == 2)
            {
                payload = double.Parse(heightparts[1]);
                if (!double.IsFinite(payload)) payload = 1;
            }

            if (double.IsFinite(b) && double.IsFinite(h))
            {
                var result = Pentate(b, h, payload);
                _sign = result._sign;
                _layer = result._layer;
                _mag = result._mag;
                return this;
            }
        }

        var tetrationparts = value.Split("^^");
        if (tetrationparts.Length == 2)
        {
            var b = double.Parse(tetrationparts[0]);
            var h = double.Parse(tetrationparts[1]);
            var heightparts = tetrationparts[1].Split(";");
            double payload = 1;

            if (heightparts.Length == 2)
            {
                payload = double.Parse(heightparts[1]);
                if (!double.IsFinite(payload)) payload = 1;
            }

            if (double.IsFinite(b) && double.IsFinite(h))
            {
                var result = Tetrate(b, h, payload);
                _sign = result._sign;
                _layer = result._layer;
                _mag = result._mag;
                return this;
            }
        }

        var powparts = value.Split("^");
        if (powparts.Length == 2)
        {
            var b = double.Parse(powparts[0]);
            var e = double.Parse(powparts[1]);

            if (double.IsFinite(b) && double.IsFinite(e))
            {
                var result = Pow(b, e);
                _sign = result._sign;
                _layer = result._layer;
                _mag = result._mag;
                return this;
            }
        }

        value = value.Trim().ToLower();
        double ba;
        double height;
        var ptparts = value.Split("pt");

        if (ptparts.Length == 2)
        {
            ba = 10;
            height = double.Parse(ptparts[0]);
            ptparts[1] = ptparts[1].Replace("(", "");
            ptparts[1] = ptparts[1].Replace(")", "");
            var payload = double.Parse(ptparts[1]);

            if (!double.IsFinite(payload)) payload = 1;
            if (double.IsFinite(ba) && double.IsFinite(height))
            {
                var result = Tetrate(ba, height, payload);
                _sign = result._sign;
                _layer = result._layer;
                _mag = result._mag;
                return this;
            }
        }

        ptparts = value.Split("p");
        if (ptparts.Length == 2)
        {
            ba = 10;
            height = double.Parse(ptparts[0]);
            ptparts[1] = ptparts[1].Replace("(", "");
            ptparts[1] = ptparts[1].Replace(")", "");
            var payload = double.Parse(ptparts[1]);

            if (!double.IsFinite(payload)) payload = 1;
            if (double.IsFinite(ba) && double.IsFinite(height))
            {
                var result = Tetrate(ba, height, payload);
                _sign = result._sign;
                _layer = result._layer;
                _mag = result._mag;
                return this;
            }
        }

        var parts = value.Split("e");
        var ecount = parts.Length - 1;
        double numberAttempt;

        switch (ecount)
        {
            case 0:
                numberAttempt = double.Parse(value);
                if (double.IsFinite(numberAttempt)) return N(numberAttempt);
                break;
            case 1:
                numberAttempt = double.Parse(value);
                if (double.IsFinite(numberAttempt) && numberAttempt != 0)
                    return N(numberAttempt);
                break;
        }

        var newparts = value.Split("e^");
        if (newparts.Length == 2)
        {
            _sign = 1;
            if (newparts[0][0] == '-') _sign = -1;
            var layerstring = "";

            for (var i = 0; i < newparts[1].Length; i++)
            {
                int chrcode = newparts[1][i];
                if ((chrcode >= 43 && chrcode <= 57) || chrcode == 101) layerstring += newparts[1][i];
                else
                {
                    _layer = long.Parse(layerstring);
                    _mag = double.Parse(newparts[1][(i + 1)..]);
                    Normalize();
                    return this;
                }
            }
        }

        if (ecount < 1)
        {
            _sign = 0;
            _layer = 0;
            _mag = 0;
            return this;
        }

        var mantissa = double.Parse(parts[0]);
        if (mantissa == 0)
        {
            _sign = 0;
            _layer = 0;
            _mag = 0;
            return this;
        }

        var exponent = double.Parse(parts[parts.Length - 1]);
        if (ecount >= 2)
        {
            var me = double.Parse(parts[parts.Length - 2]);
            if (double.IsFinite(me))
            {
                exponent *= Math.Sign(me);
                exponent += FMagLog10(me);
            }
        }

        if (!double.IsFinite(mantissa))
        {
            _sign = parts[0] == "-" ? -1 : 1;
            _layer = ecount;
            _mag = exponent;
        }
        else if (ecount == 1)
        {
            _sign = Math.Sign(mantissa);
            _layer = 1;
            _mag = exponent + Math.Log10(Math.Abs(mantissa));
        }
        else
        {
            _sign = Math.Sign(mantissa);
            _layer = ecount;
            if (ecount == 2)
            {
                var result = Mul(FC(1, 2, exponent), D(mantissa));
                _sign = result._sign;
                _layer = result._layer;
                _mag = result._mag;
                return this;
            }

            _mag = exponent;
        }

        Normalize();
        return this;
    }

    public Eternity Abs() => FC_NN(_sign == 0 ? 0 : 1, _layer, _mag);
    public Eternity Neg() => FC_NN(-_sign, _layer, _mag);

    public Eternity Round()
    {
        return _mag < 0
            ? Zero
            : _layer == 0
                ? FC(_sign, 0, Math.Round(_mag))
                : this;
    }

    public Eternity Floor()
    {
        return _mag < 0
            ? Zero
            : _layer == 0
                ? FC(_sign, 0, Math.Floor(_mag))
                : this;
    }

    public Eternity Ceil()
    {
        return _mag < 0
            ? Zero
            : _layer == 0
                ? FC(_sign, 0, Math.Ceiling(_mag))
                : this;
    }

    public Eternity Trunc()
    {
        return _mag < 0
            ? Zero
            : _layer == 0
                ? FC(_sign, 0, Math.Truncate(_mag))
                : this;
    }

    public Eternity Add(Eternity d)
    {
        if (!double.IsFinite(_layer)) return this;
        if (!double.IsFinite(d._layer)) return d;
        if (_sign == 0) return d;
        if (d._sign == 0) return this;
        if (_sign == -d._sign && _layer == d._layer && _mag == d._mag) return Zero;
        if (_layer >= 2 || d._layer >= 2) return AbsMax(d);

        Eternity a, b;
        if (d.CmpAbs(this) > 0)
        {
            a = this;
            b = d;
        }
        else
        {
            a = d;
            b = this;
        }

        if (a._layer == 0 && b._layer == 0) return DE(a._sign * a._mag + b._sign * b._mag);
        double layera = a._layer * Math.Sign(a._mag);
        double layerb = b._layer * Math.Sign(b._mag);

        if (layera - layerb >= 2) return a;

        double mantissa;
        double magdiff;
        switch (layera)
        {
            case 0 when layerb == -1:
                if (Math.Abs(b._mag - Math.Log10(a._mag)) > MaxSignificantDigits) return a;
                magdiff = Math.Pow(10, Math.Log10(a._mag) - b._mag);
                mantissa = b._sign + a._sign * magdiff;
                return FC(Math.Sign(mantissa), 1, b._mag + Math.Log10(Math.Abs(mantissa)));
            case 1 when layerb == 0:
                if (Math.Abs(a._mag - Math.Log10(b._mag)) > MaxSignificantDigits) return a;
                magdiff = Math.Pow(10, a._mag - Math.Log10(b._mag));
                mantissa = b._sign + a._sign * magdiff;
                return FC(Math.Sign(mantissa), 1, Math.Log10(b._mag) + Math.Log10(Math.Abs(mantissa)));
        }

        if (Math.Abs(a._mag - b._mag) > MaxSignificantDigits) return a;
        magdiff = Math.Pow(10, a._mag - b._mag);
        mantissa = b._sign + a._sign * magdiff;
        return FC(Math.Sign(mantissa), 1, b._mag + Math.Log10(Math.Abs(mantissa)));
    }

    public Eternity Sub(Eternity d) => Add(D(d).Neg());

    public Eternity Mul(Eternity d)
    {
        if (!double.IsFinite(_mag)) return D(this);
        if (!double.IsFinite(d._mag)) return D(d);
        if (_sign == 0 || d._sign == 0) return FC_NN(0, 0, 0);
        if (_layer == d._layer && _mag == -d._mag) return FC_NN(_sign * d._sign, 0, 1);

        Eternity a, b;
        if (_layer > d._layer || (_layer == d._layer && Math.Abs(_mag) > Math.Abs(d._mag)))
        {
            a = this;
            b = d;
        }
        else
        {
            a = d;
            b = this;
        }

        if (a._layer == 0 && b._layer == 0) return DE(a._sign * b._sign * a._mag * b._mag);
        if (a._layer >= 3 || a._layer - b._layer >= 2) return FC(a._sign * b._sign, a._layer, a._mag);
        Eternity newmag;

        switch (a._layer)
        {
            case 1 when b._layer == 0:
                return FC(a._sign * b._sign, 1, a._mag + Math.Log10(b._mag));
            case 1 when b._layer == 1:
                return FC(a._sign * b._sign, 1, a._mag + b._mag);
            case 2 when b._layer == 1:
                newmag = FC(Math.Sign(a._mag), a._layer - 1, Math.Abs(a._mag))
                    .Add(FC(Math.Sign(b._mag), b._layer - 1, Math.Abs(b._mag)));
                return FC(a._sign * b._sign, newmag._layer + 1, newmag._sign * newmag._mag);
            case 2 when b._layer == 2:
                newmag = FC(Math.Sign(a._mag), a._layer - 1, Math.Abs(a._mag))
                    .Add(FC(Math.Sign(b._mag), b._layer - 1, Math.Abs(b._mag)));
                return FC(a._sign * b._sign, newmag._layer + 1, newmag._sign * newmag._mag);
            default:
                return NaN;
        }
    }

    public Eternity Div(Eternity d) => Mul(d.Recip());

    public Eternity Recip()
    {
        return _mag == 0
            ? NaN
            : _layer == 0
                ? FC(_sign, 0, 1 / _mag)
                : FC(_sign, _layer, -_mag);
    }

    public Eternity Max(Eternity d) => Lt(d) ? d : this;
    public Eternity Min(Eternity d) => Gt(d) ? d : this;
    public Eternity AbsMax(Eternity d) => CmpAbs(d) < 0 ? d : this;
    public Eternity AbsMin(Eternity d) => CmpAbs(d) > 0 ? d : this;
    public Eternity Clamp(Eternity d1, Eternity d2) => Max(d1).Min(d2);

    public Eternity AbsLog10()
    {
        return _sign == 0
            ? NaN
            : _layer > 0
                ? FC(Math.Sign(_mag), _layer - 1, Math.Abs(_mag))
                : FC(1, 0, Math.Log10(_mag));
    }

    public Eternity Log10()
    {
        return _sign <= 0
            ? NaN
            : _layer > 0
                ? FC(Math.Sign(_mag), _layer - 1, Math.Abs(_mag))
                : FC(_sign, 0, Math.Log10(_mag));
    }

    public Eternity Log() => Log(MathE);

    public Eternity Log(Eternity d)
    {
        return _sign <= 0
            ? NaN
            : d._sign <= 0
                ? NaN
                : d._sign == 1 && d._layer == 0 && d._mag == 1
                    ? NaN
                    : _layer == 0 && d._layer == 0
                        ? FC(_sign, 0, Math.Log(_mag) / Math.Log(d._mag))
                        : Div(Log10(), d.Log10());
    }

    public Eternity Log2()
    {
        return _sign <= 0
            ? NaN
            : _layer == 0
                ? FC(_sign, 0, Math.Log2(_mag))
                : _layer == 1
                    ? FC(Math.Sign(_mag), 0, Math.Abs(_mag) * 3.321928094887362)
                    : _layer == 2
                        ? FC(Math.Sign(_mag), 1, Math.Abs(_mag) + 0.5213902276543247)
                        : FC(Math.Sign(_mag), _layer - 1, Math.Abs(_mag));
    }

    public Eternity Ln()
    {
        return _sign <= 0
            ? NaN
            : _layer == 0
                ? FC(_sign, 0, Math.Log(_mag))
                : _layer == 1
                    ? FC(Math.Sign(_mag), 0, Math.Abs(_mag) * 2.302585092994046)
                    : _layer == 2
                        ? FC(Math.Sign(_mag), 1, Math.Abs(_mag) + 0.36221568869946325)
                        : FC(Math.Sign(_mag), _layer - 1, Math.Abs(_mag));
    }

    public Eternity Pow(Eternity d)
    {
        var a = this;
        var b = d;
        switch (a._sign)
        {
            case 0:
                return b.Eq(0) ? FC_NN(1, 0, 1) : a;
            case 1 when a._layer == 0 && a._mag == 1:
                return a;
        }

        switch (b._sign)
        {
            case 0:
                return FC_NN(1, 0, 1);
            case 1 when b._layer == 0 && b._mag == 1:
                return a;
        }

        var result = a.AbsLog10().Mul(b).Pow10();
        if (_sign != -1) return result;

        if (Math.Abs(b.ToNumber() % 2) % 2 == 1) return result.Neg();
        return Math.Abs(b.ToNumber() % 2) % 2 == 0 ? result : NaN;
    }

    public Eternity Pow10()
    {
        if (!double.IsFinite(_layer) || !double.IsFinite(_mag)) return NaN;
        var a = this;
        if (a._layer == 0)
        {
            var newmag = Math.Pow(10, a._sign * a._mag);
            if (double.IsFinite(newmag) && Math.Abs(newmag) >= 0.1) return FC(1, 0, newmag);
            if (a._sign == 0) return One;
            a = FC_NN(a._sign, a._layer + 1, Math.Log10(a._mag));
        }

        return a._sign switch
        {
            > 0 when a._mag >= 0 => FC(a._sign, a._layer + 1, a._mag),
            < 0 when a._mag >= 0 => FC(-a._sign, a._layer + 1, -a._mag),
            _ => One
        };
    }

    public Eternity Root(Eternity d) => Pow(d.Recip());

    public Eternity Factorial()
    {
        if (_mag < 0) return Add(1).Gamma();
        return _layer switch
        {
            0 => Add(1).Gamma(),
            1 => Exp(Mul(this, Ln(this).Sub(1))),
            _ => Exp(this)
        };
    }

    public Eternity Gamma()
    {
        if (_mag < 0) return Recip();

        switch (_layer)
        {
            case 0 when Lt(FC_NN(1, 0, 24)):
                return N(FGamma(_sign * _mag));
            case 0:
            {
                var t = _mag - 1;
                var l = 0.9189385332046727;
                l += (t + 0.5) * Math.Log(t);
                l -= t;
                var n2 = t * t;
                var np = t;
                var lm = 12 * np;
                var adj = 1 / lm;
                var l2 = l + adj;
                if (l2 == l) return Exp(l);

                l = l2;
                np *= n2;
                lm = 360 * np;
                adj = 1 / lm;
                l2 = l - adj;
                if (l2 == l) return Exp(l);

                l = l2;
                np *= n2;
                lm = 1260 * np;
                var lt = 1 / lm;
                l += lt;
                np *= n2;
                lm = 1680 * np;
                lt = 1 / lm;
                l -= lt;
                return Exp(l);
            }
            case 1:
                return Exp(Mul(this, Ln(this).Sub(1)));
            default:
                return Exp(this);
        }
    }

    public Eternity LnGamma() => Gamma().Ln();

    public Eternity Exp()
    {
        if (_mag < 0) return One;

        return _layer switch
        {
            0 when _mag <= 709.7 => N(Math.Exp(_sign * _mag)),
            0 => FC(1, 1, _sign * Math.Log10(Math.E) * _mag),
            1 => FC(1, 2, _sign * (Math.Log10(0.4342944819032518) + _mag)),
            _ => FC(1, _layer + 1, _sign * _mag)
        };
    }

    public Eternity Sq() => Pow(2);

    public Eternity Sqrt()
    {
        switch (_layer)
        {
            case 0:
                return N(Math.Sqrt(_sign * _mag));
            case 1:
                return FC(1, 2, Math.Log10(_mag) - 0.3010299956639812);
        }

        var result = Div(FC_NN(_sign, _layer - 1, _mag), FC_NN(1, 0, 2));
        result._layer += 1;
        result.Normalize();
        return result;
    }

    public Eternity Cube() => Pow(3);
    public Eternity Cbrt() => Pow(1 / 3);

    public Eternity Tetrate(double height = 2, Eternity? payload = null)
    {
        payload ??= FC_NN(1, 0, 1);
        switch (height)
        {
            case 1:
                return Pow(this, payload);
            case 0:
                return D(payload);
        }

        if (Eq(One)) return One;
        if (Eq(-1)) return Pow(this, payload);
        if (height == double.PositiveInfinity)
        {
            var this_num = ToNumber();

            switch (this_num)
            {
                case <= 1.44466786100976613366 and >= 0.06598803584531253708 when this_num > 1.444667861009099:
                    return N(Math.E);
                case <= 1.44466786100976613366 and >= 0.06598803584531253708:
                    var negln = Ln(this).Neg();
                    return negln.LambertW().Div(negln);
                case > 1.44466786100976613366:
                    return N(double.PositiveInfinity);
                default:
                    return NaN;
            }
        }

        if (Eq(Zero))
        {
            var result = Math.Abs((height + 1) % 2);
            if (result > 1) result = 2 - result;
            return N(result);
        }

        if (height < 0) return IteratedLog(payload, this, -height);

        payload = D(payload);
        var oldheight = height;
        height = Math.Truncate(height);
        var fracheight = oldheight - height;

        if (Gt(Zero) && Lte(1.44466786100976613366))
        {
            height = Math.Min(10000, height);
            for (var i = 0; i < height; i++)
            {
                var old_payload = payload;
                payload = Pow(payload);
                if (old_payload.Eq(payload)) return payload;
            }

            if (fracheight == 0) return payload;
            var next_payload = Pow(payload);
            return payload.Mul(1 - fracheight).Add(next_payload.Mul(fracheight));
        }

        if (fracheight != 0)
        {
            if (payload.Eq(One))
            {
                if (Gt(10)) payload = Pow(fracheight);
                else
                {
                    payload = N(TetrateCritical(ToNumber(), fracheight));
                    if (Lt(2)) payload = payload.Sub(1).Mul(Sub(1)).Add(1);
                }
            }
            else
            {
                if (Eq(10)) payload = payload.LayerAdd10(fracheight);
                else payload = payload.LayerAdd(fracheight, this);
            }
        }

        for (var i = 0; i < height; i++)
        {
            payload = Pow(payload);
            if (!double.IsFinite(payload._mag)) return payload.Normalize();
            if (payload._layer - _layer > 3)
            {
                return FC_NN(payload._sign, (long) (payload._layer + (height - i - 1)), payload._mag);
            }

            if (i > 10000) return payload;
        }

        return payload;
    }

    public Eternity IteratedLog(Eternity? b = null, double times = 1)
    {
        b ??= 10;
        if (times < 0) return Tetrate(b, -times, this);
        var result = D(this);
        var fulltimes = times;
        times = Math.Truncate(times);
        var fraction = fulltimes - times;

        if (result._layer - b._layer > 3)
        {
            var layerloss = (long) Math.Min(times, result._layer - b._layer - 3);
            times -= layerloss;
            result._layer -= layerloss;
        }

        for (var i = 0; i < times; i++)
        {
            result = result.Log(b);
            if (!double.IsFinite(result._mag)) return result.Normalize();
            if (i > 10000) return result;
        }

        if (fraction is > 0 and < 1)
        {
            result = b.Eq(10) ? result.LayerAdd10(-fraction) : result.LayerAdd(-fraction, b);
        }

        return result;
    }

    public Eternity Slog(Eternity? b)
    {
        b ??= 10;
        if (b.Lte(Zero)) return NaN;
        if (b.Eq(One)) return NaN;
        if (b.Lt(One))
        {
            if (Eq(One)) return Zero;
            return Eq(Zero) ? NegativeOne : NaN;
        }

        if (_mag < 0 || Eq(Zero)) return NegativeOne;
        double result = 0;
        var copy = D(this);
        if (copy._layer - b._layer > 3)
        {
            var layerloss = copy._layer - b._layer - 3;
            result += layerloss;
            copy._layer -= layerloss;
        }

        for (var i = 0; i < 100; i++)
        {
            if (copy.Lt(Zero))
            {
                copy = Pow(b, copy);
                result -= 1;
            }
            else if (copy.Lte(One)) return N(result + SLogCritical(b.ToNumber(), copy.ToNumber()));
            else
            {
                result += 1;
                copy = Log(copy, b);
            }
        }

        return N(result);
    }

    public Eternity LayerAdd10(double diff)
    {
        var result = D(this);
        if (diff >= 1)
        {
            if (result._mag < 0 && result._layer > 0)
            {
                result._sign = 0;
                result._mag = 0;
                result._layer = 0;
            }
            else if (result._sign == -1 && result._layer == 0)
            {
                result._sign = 1;
                result._mag = -result._mag;
            }

            var layeradd = (long) Math.Truncate(diff);
            diff -= layeradd;
            result._layer += layeradd;
        }

        if (diff <= -1)
        {
            var layeradd = (long) Math.Truncate(diff);
            diff -= layeradd;
            result._layer += layeradd;
            if (result._layer < 0)
            {
                for (var i = 0; i < 100; i++)
                {
                    result._layer++;
                    result._mag = Math.Log10(result._mag);
                    if (!double.IsFinite(result._mag))
                    {
                        if (result._sign == 0) result._sign = 1;
                        if (result._layer < 0) result._layer = 0;
                        return result.Normalize();
                    }

                    if (result._layer >= 0) break;
                }
            }
        }

        while (result._layer < 0)
        {
            result._layer++;
            result._mag = Math.Log10(result._mag);
        }

        if (result._sign == 0)
        {
            result._sign = 1;
            if (result._mag == 0 && result._layer >= 1)
            {
                result._layer -= 1;
                result._mag = 1;
            }
        }

        result.Normalize();
        return diff != 0 ? result.LayerAdd(diff, 10) : result;
    }

    public Eternity LayerAdd(double diff, Eternity b)
    {
        var slogthis = SLog(b).ToNumber();
        var slogdest = slogthis + diff;
        if (slogdest >= 0) return Tetrate(b, slogdest);

        return !double.IsFinite(slogdest)
            ? NaN
            : Log(slogdest >= -1 ? Tetrate(b, slogdest + 1) : Log(Tetrate(b, slogdest + 2), b), b);
    }

    public Eternity LambertW()
    {
        if (Lt(-0.3678794411710499)) throw new Exception("lambertw is unimplemented for results less than -1, sorry!");
        if (_mag < 0) return N(FLambertW(ToNumber()));

        return _layer switch
        {
            0 => N(FLambertW(_sign * _mag)),
            1 => DLambertW(this),
            2 => DLambertW(this),
            >= 3 => FC_NN(_sign, _layer - 1, _mag),
            _ => throw new Exception("Unhandled behavior in lambertw()")
        };
    }

    public Eternity SSqrt()
    {
        if (_sign == 1 && _layer >= 3) return FC_NN(_sign, _layer - 1, _mag);
        var lnx = Ln();
        return lnx.Div(lnx.LambertW());
    }

    public Eternity Pentate(double height = 2, Eternity? payload = null)
    {
        payload ??= FC_NN(1, 0, 1);
        var oldheight = height;
        height = Math.Truncate(height);
        var fracheight = oldheight - height;
        if (fracheight != 0)
        {
            if (payload.Eq(One))
            {
                height++;
                payload = N(fracheight);
            }
            else payload = Eq(10) ? payload.LayerAdd10(fracheight) : payload.LayerAdd(fracheight, this);
        }

        for (var i = 0; i < height; i++)
        {
            payload = Tetrate(payload.ToNumber());
            if (!double.IsFinite(payload._mag)) return payload.Normalize();
            if (i > 10) return payload;
        }

        return payload;
    }

    public Eternity Sin()
    {
        return _mag < 0
            ? this
            : _layer == 0
                ? N(Math.Sin(_sign * _mag))
                : FC_NN(0, 0, 0);
    }

    public Eternity Cos()
    {
        return _mag < 0
            ? One
            : _layer == 0
                ? N(Math.Cos(_sign * _mag))
                : FC_NN(0, 0, 0);
    }

    public Eternity Tan()
    {
        return _mag < 0
            ? this
            : _layer == 0
                ? N(Math.Tan(_sign * _mag))
                : FC_NN(0, 0, 0);
    }

    public Eternity Asin()
    {
        return _mag < 0
            ? this
            : _layer == 0
                ? N(Math.Asin(_sign * _mag))
                : FC_NN(0, 0, double.NaN);
    }

    public Eternity Acos()
    {
        return _mag < 0
            ? N(Math.Acos(ToNumber()))
            : _layer == 0
                ? N(Math.Acos(_sign * _mag))
                : FC_NN(0, 0, double.NaN);
    }

    public Eternity Atan()
    {
        return _mag < 0
            ? this
            : _layer == 0
                ? N(Math.Atan(_sign * _mag))
                : N(Math.Atan(_sign * 1.79e308));
    }

    public Eternity Sinh() => Exp().Sub(Neg().Exp()).Div(2);
    public Eternity Cosh() => Exp().Add(Neg().Exp()).Div(2);
    public Eternity Tanh() => Sinh().Div(Cosh());
    public Eternity Asinh() => Ln(Add(Sq().Add(1).Sqrt()));
    public Eternity Acosh() => Ln(Add(Sq().Sub(1).Sqrt()));
    public Eternity Atanh() => Abs().Gte(1) ? FC_NN(0, 0, double.NaN) : Ln(Add(1).Div(N(1).Sub(this))).Div(2);

    public int Sign() => _sign;

    public int Cmp(Eternity d)
    {
        return _sign > d._sign
            ? 1
            : _sign < d._sign
                ? -1
                : _sign * CmpAbs(d);
    }

    public int CmpAbs(Eternity d)
    {
        double layera = _mag > 0 ? _layer : -_layer;
        double layerb = d._mag > 0 ? d._layer : -d._layer;
        if (layera > layerb) return 1;
        if (layera < layerb) return -1;

        if (_mag > d._mag) return 1;
        if (_mag < d._mag) return -1;
        return 0;
    }

    public bool IsNan() => double.IsNaN(_mag);
    public bool IsFinite() => double.IsFinite(_mag);
    public bool Eq(Eternity d) => _sign == d._sign && _layer == d._layer && _mag == d._mag;
    public bool Neq(Eternity d) => !Eq(d);
    public bool Lt(Eternity d) => Cmp(d) == -1;
    public bool Lte(Eternity d) => !Gt(d);
    public bool Gt(Eternity d) => Cmp(d) == 1;
    public bool Gte(Eternity d) => !Lt(d);

    public double ToNumber()
    {
        if (!double.IsFinite(_layer)) return double.NaN;

        return _layer switch
        {
            0 => _sign * _mag,
            1 => _sign * Math.Pow(10, _mag),
            _ => _mag > 0 ? _sign > 0 ? double.PositiveInfinity : double.NegativeInfinity : 0
        };
    }

    public double MantissaWithDecPlaces(int places)
    {
        return M switch
        {
            double.NaN => double.NaN,
            0 => 0,
            _ => DecPlaces(M, places)
        };
    }

    public double MagnitudeWithDecPlaces(int places)
    {
        return _mag switch
        {
            double.NaN => double.NaN,
            0 => 0,
            _ => DecPlaces(_mag, places)
        };
    }

    public override string ToString()
    {
        if (double.IsNaN(_mag)) return "NaN";
        if (_mag == double.PositiveInfinity || _layer == long.MinValue) return _sign == 1 ? "Infinity" : "-Infinity";

        return _layer switch
        {
            0 when (_mag < 1e21 && _mag > 1e-7) || _mag == 0 => (_sign * _mag).ToString(CultureInfo.InvariantCulture),
            0 => $"{M}e{E}",
            1 => $"{M}e{E}",
            _ => _layer <= MaxEsInARow
                ? $"{(_sign == -1 ? "-" : "")}{new string('e', (int) _layer)}{_mag}"
                : $"{(_sign == -1 ? "-" : "")}(e^{_layer}){_mag}"
        };
    }

    public string ToExponential(int places)
    {
        return _layer == 0 ? N(_sign * _mag).ToExponential(places) : ToStringWithDecPlaces(places);
    }

    public string ToFixed(int places)
    {
        return _layer == 0 ? N(_sign * _mag).ToFixed(places) : ToStringWithDecPlaces(places);
    }

    public string ToPrecision(int places)
    {
        if (E <= -7) return ToExponential(places - 1);
        return places > E ? ToFixed((int) (places - E - 1)) : ToExponential(places - 1);
    }

    public string ToStringWithDecPlaces(int places)
    {
        return _layer switch
        {
            0 when (_mag < 1e21 && _mag > 1e-7) || _mag == 0 => N(_sign * _mag).ToFixed(places),
            0 => $"{DecPlaces(M, places)}e{DecPlaces(E, places)}",
            1 => $"{DecPlaces(M, places)}e{DecPlaces(E, places)}",
            _ => _layer <= MaxEsInARow
                ? $"{(_sign == -1 ? "-" : "")}{new string('e', (int) _layer)}{DecPlaces(_mag, places)}"
                : $"{(_sign == -1 ? "-" : "")}(e^{_layer}){DecPlaces(_mag, places)}"
        };
    }

    public string FormatDec(int? p = null, Format? f = null)
    {
        f ??= DecFormat;
        return f switch
        {
            Format.Mixed => FormatStandard(p ?? 2),
            _ => "Invalid Format",
        };
    }

    public string FormatStandard(int p = 2)
    {
        if (Lt(Abs(), 1e3))
        {
            return
                $"{(CmpAbs(Zero) < 0 ? "-" : "")}{Math.Abs(_mag).ToString($"##0{(p > 0 ? $".{new string('0', p)}" : "")}")}";
        }

        return Lt(Abs(), 1e9)
            ? $"{(CmpAbs(Zero) < 0 ? "-" : "")}{Math.Abs(_mag):###,###,###}"
            : $"{(CmpAbs(Zero) < 0 ? "-" : "")}{Math.Abs(M).ToString($"0{(p > 0 ? $".{new string('0', p)}" : "")}")}e{Math.Abs(E):###,###,##0}";
    }

    public override bool Equals(object? obj) => ReferenceEquals(this, obj);
    public override int GetHashCode() => throw new NotImplementedException();
}