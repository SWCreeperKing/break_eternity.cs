// See https://aka.ms/new-console-template for more information

using Break_Eternity;
using NUnit.Framework;

public class Tester
{
    [TestFixture]
    public class TestClass
    {
        public Eternity n1;
        public Eternity n2;

        [SetUp]
        public void Setup()
        {
            n1 = 10000;
            // n2 = "1e50";
            n2 = new Eternity(1, 50);
            Console.WriteLine(n2);
        }

        [Test]
        public void NullTest()
        {
            Eternity nc1 = 6;
            Eternity? nc2 = null;
            Assert.IsNotNull(nc1);
            Assert.IsNull(nc2);
        }

        [Test]
        public void PlusPlus()
        {
            Eternity n = 1;
            Assert.True(n == 1);
            n++;
            Assert.True(n == 2 && ++n == 3);
        }

        [Test]
        public void LogTest()
        {
            Assert.True(n1.Log10() == 4 && n2.Log10() == 50 && n2.Log(n1) == 12.5);
            Assert.True((double)n2.Log() == 115.12925464970229);
        }

        [Test]
        public void AddTest()
        {
            Console.WriteLine(n1 + n2);
            Console.WriteLine(new Eternity(1, 50));
            Assert.True(n1 + n2 == new Eternity(1, 50));
            var prepreadd = new Eternity("1e30");
            var preAdd = prepreadd + "6e29";
            var addRes = new Eternity(1.6f, 30);
            Assert.True(preAdd == addRes);
        }

        [Test]
        public void SubTest()
        {
            var sub1 = n1 - n2;
            var sub1Res = new Eternity(-1, 50);
            var sub2 = n2 - n1;
            var sub2Res = new Eternity(1, 50);
            Assert.True(sub1 == sub1Res);
            Assert.True(sub2 == sub2Res);
            Assert.True(Eternity.One / 2 == .5f);
            var sub3n1 = new Eternity(2.2561784573198107, 6);
            var sub3n2 = new Eternity(1.9603208295529388, 7);
            var sub3 = sub3n1 - sub3n2;
            var sub3Res = new Eternity(-1.73470298382095773, 7);
            Assert.True(sub3 == sub3Res);
        }

        [Test]
        public void MultiTest()
        {
            Assert.True(n1 * n2 == new Eternity(1, 54));
        }

        [Test]
        public void DivTest()
        {
            Assert.True(n1 / n2 == new Eternity(1, -46));
            Assert.True(n2 / n1 == new Eternity(1, 46));
        }

        [Test]
        public void PowerTest()
        {
            Assert.True(n1.Pow(n2) == new Eternity(1, 4e50));
            Assert.True(n2.Pow(n1) == new Eternity(1, 50e4));
        }

        [Test]
        public void CompareTest()
        {
            Assert.False(n1 > n2);
            Assert.True(n1 < n2);
            Assert.True(n1 != n2);
            Assert.False(n1 == n2);
        }

        [Test]
        public void RootTest()
        {
            Assert.True(new Eternity(3).Pow(4).Root(4) == 3);
        }
        
        /*
        Dec b = 2;

        var n = new Dec(1e12) - new Dec(9.2e11);
        Console.WriteLine(n);

        // Dec.formatter = Engineering.Name;
        SciNotation.beforeSciCut = 8;
        Engineering.beforeSciCut = 8;

        var convert = JsonConvert.SerializeObject(new Dec(5.342, 3));
        Console.WriteLine("e");
        Console.WriteLine(convert);
        Console.WriteLine(JsonConvert.DeserializeObject<Dec>(convert));
        Console.WriteLine("E");
         */
    }
} 