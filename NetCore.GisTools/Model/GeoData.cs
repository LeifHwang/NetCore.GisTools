using System.Diagnostics;

namespace NetCore.GisTools.Model
{
    [DebuggerDisplay("({Longitude},{Latitude}) {Value}")]
    public class GeoData<T> : GeoPoint
    {
        public T Value { get; set; }

        public GeoData(double x, double y, T value) : base(x, y)
        {
            Value = value;
        }
    }

    public class GeoDoubleData : GeoData<double>
    {
        public GeoDoubleData(double x, double y, double value) : base(x, y, value) { }
    }
    public class GeoDecimalData : GeoData<decimal>
    {
        public GeoDecimalData(double x, double y, decimal value) : base(x, y, value) { }
    }
}
