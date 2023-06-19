using System.Collections.Generic;
using System.Diagnostics;

namespace NetCore.GisTools.Demo.Model
{
    public class KrigingInput
    {
        public IEnumerable<PointData> points { get; set; }
        public int sigma2 { get; set; }
        public int alpha { get; set; }

        public IEnumerable<Polygon> polygons { get; set; }
        public double width { get; set; }
    }

    [DebuggerDisplay("({x},{y})")]
    public class Point
    {
        public double x { get; set; }
        public double y { get; set; }
    }
    [DebuggerDisplay("({x},{y}) {value}")]
    public class PointData : Point
    {
        public double value { get; set; }
    }

    public class Polygon
    {
        public IEnumerable<Point> points { get; set; }
    }
}
