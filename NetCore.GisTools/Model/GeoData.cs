using System.Diagnostics;

namespace NetCore.GisTools.Model
{
    [DebuggerDisplay("({Longitude},{Latitude}) {Value}")]
    public class GeoData<T>
    {
        public T Value { get; set; }

        public double Longitude { get; set; }
        public double Latitude { get; set; }
    }

    public class GeoDoubleData : GeoData<double> { }
    public class GeoDecimalData : GeoData<decimal> { }
}
