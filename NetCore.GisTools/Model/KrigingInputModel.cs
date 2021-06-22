using System;
using System.Collections.Generic;
using System.Text;

namespace NetCore.GisTools.Model
{
    public class KrigingInputModel
    {
        public double[] T { get; set; }
        public double[] X { get; set; }
        public double[] Y { get; set; }

        public double MinX { get; set; }
        public double MaxX { get; set; }

        public double MinY { get; set; }
        public double MaxY { get; set; }
    }
}
