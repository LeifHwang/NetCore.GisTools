using NetCore.GisTools.Model;
using NetCore.GisTools.Utils;

using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Dynamic;
using System.Linq;

namespace NetCore.GisTools
{
    public enum KrigingTrainModel
    {
        Gaussian,
        Exponential,
        Spherical
    }

    public static class Kriging
    {
        /// <summary>
        /// Train using gaussian processes with bayesian priors
        /// </summary>
        /// <param name="values"></param>
        /// <param name="sigma2">the variance parameters of the gaussian process</param>
        /// <param name="alpha">the prior of the variogram model</param>
        /// <param name="model"></param>
        /// <returns></returns>
        public static KrigingTrainVariogram Train(GeoDoubleData[] values, int sigma2, int alpha,
            KrigingTrainModel model = KrigingTrainModel.Exponential)
        {
            if (values == null || values.Length == 0)
                throw new ArgumentNullException(nameof(values));
            if (alpha <= 0)
                throw new ArgumentException("alpha must larger then 0");

            // Lag distance/semivariance
            var inputPointCount = values.Length;
            var inputValues = new double[inputPointCount];
            double valueMin = double.MaxValue, valueMax = double.MinValue;

            var arithmeticSum = (inputPointCount * inputPointCount - inputPointCount) / 2;
            var distance = new double[arithmeticSum][];
            for (int i = 0, k = 0; i < inputPointCount; i++)
            {
                inputValues[i] = values[i].Value;
                valueMin = Math.Min(valueMin, values[i].Value);
                valueMax = Math.Max(valueMax, values[i].Value);

                for (int j = 0; j < i; j++, k++)
                {
                    var temp = new double[2];
                    temp[0] = values[i].Distance(values[j]);
                    temp[1] = Math.Abs(values[i].Value - values[j].Value);

                    distance[k] = temp;
                }
            }
            Array.Sort(distance, (a, b) => a[0].CompareTo(b[0]));

            var result = new KrigingTrainVariogram(values, model) { MinValue = valueMin, MaxValue = valueMax };

            // Bin lag distance
            var lagsCount = Math.Min(30, arithmeticSum);
            var lags = new double[lagsCount];
            var semis = new double[lagsCount];
            var binCount = 0;
            if (arithmeticSum < 30)
            {
                for (; binCount < arithmeticSum; binCount++)
                {
                    lags[binCount] = distance[binCount][0];
                    semis[binCount] = distance[binCount][1];
                }
            }
            else
            {
                var tolerance = distance.Last()[0] / 30;

                for (int i = 0, j = 0, k = 0; i < 30 && j < arithmeticSum; i++, k = 0)
                {
                    while (distance[j][0] <= ((i + 1) * tolerance))
                    {
                        lags[binCount] += distance[j][0];
                        semis[binCount] += distance[j][1];
                        j++; k++;

                        if (j >= arithmeticSum) break;
                    }

                    if (k > 0)
                    {
                        lags[binCount] /= k;
                        semis[binCount] /= k;
                        binCount++;
                    }
                }
                if (binCount < 2) return result; // Error: Not enough points
            }
            result.Range = lags[binCount - 1] - lags[0];

            // Feature transformation
            var lagModels = new double[binCount * 2];
            var y = new double[binCount];
            for (var i = 0; i < binCount; i++)
            {
                lagModels[i * 2] = 1.0;
                lagModels[i * 2 + 1] = result.Transform(lags[i]);

                y[i] = semis[i];
            }

            // Least squares
            var lagModelsTrans = Matrix.Transpose(lagModels, binCount, 2);
            var lagModelRanges = Matrix.Multiply(lagModelsTrans, lagModels, 2, binCount, 2);
            lagModelRanges = Matrix.Add(lagModelRanges, Matrix.Diag(1 / alpha, 2), 2);

            var cloneLaglagModelRanges = Matrix.Clone(lagModelRanges);

            if (Matrix.Chol(lagModelRanges, 2))
                Matrix.Chol2Inv(lagModelRanges, 2);
            else
            {
                Matrix.Solve(cloneLaglagModelRanges, 2);
                lagModelRanges = cloneLaglagModelRanges;
            }
            var rectangle = Matrix.Multiply(Matrix.Multiply(lagModelRanges, lagModelsTrans, 2, 2, binCount)
                , y, 2, binCount, 1);

            // Variogram parameters
            result.Nugget = rectangle[0];
            result.Sill = rectangle[1] * result.Range + rectangle[0];

            // Gram matrix with prior
            var models = new double[inputPointCount * inputPointCount];
            for (var i = 0; i < inputPointCount; i++)
            {
                for (var j = 0; j < i; j++)
                {
                    var modelItem = result.Model(values[i].Distance(values[j]));

                    models[i * inputPointCount + j] = modelItem;
                    models[j * inputPointCount + i] = modelItem;
                }
                models[i * inputPointCount + i] = result.Model(0);
            }

            // Inverse penalized Gram matrix projected to target vector
            var sigmaModels = Matrix.Add(models, Matrix.Diag(sigma2, inputPointCount), inputPointCount);
            var cloneSigmaModels = Matrix.Clone(sigmaModels);

            if (Matrix.Chol(sigmaModels, inputPointCount))
                Matrix.Chol2Inv(sigmaModels, inputPointCount);
            else
            {
                Matrix.Solve(cloneSigmaModels, inputPointCount);
                sigmaModels = cloneSigmaModels;
            }

            // Copy unprojected inverted matrix as K
            result.K = Matrix.Clone(sigmaModels);
            result.M = Matrix.Multiply(sigmaModels, inputValues, inputPointCount, inputPointCount, 1);

            return result;
        }

        /// <summary>
        /// Gridded matrices or contour paths
        /// <param name="width">the grid size</param>
        /// </summary>
        public static KigingGrid Grid(GeoPolygon[] polygons, KrigingTrainVariogram variogram, double width)
        {
            if (width <= 0)
                throw new ArgumentException("width must larger then 0");

            if (polygons.Length == 0)
                return new KigingGrid();

            // Boundaries of polygons space
            var maxOuterExtend = polygons[0].GetExtend();
            for (int i = 1; i < polygons.Length; i++)
                maxOuterExtend.Union(polygons[i].GetExtend());

            // Alloc for O(n^2) space
            var grid = new double?[(int)Math.Ceiling(maxOuterExtend.Width / width) + 1][];
            for (var i = 0; i < grid.Length; i++)
                grid[i] = new double?[(int)Math.Ceiling(maxOuterExtend.Height / width) + 1];

            var values = new List<KigingGridValue>();
            foreach (var polygon in polygons)
            {
                var extend = polygon.GetExtend();

                var rowMin = (int)Math.Floor(((extend.LngMin - maxOuterExtend.LngMin) * (1 - 1 % width) / width));
                var rowMax = (int)Math.Ceiling((extend.LngMax - ((extend.LngMax - maxOuterExtend.LngMax) % width) -
                    maxOuterExtend.LngMin) / width);

                var colMin = (int)Math.Floor(((extend.LatMin - maxOuterExtend.LatMin) * (1 - 1 % width) / width));
                var colMax = (int)Math.Ceiling((extend.LatMax - ((extend.LatMax - maxOuterExtend.LatMax) % width) -
                    maxOuterExtend.LatMin) / width);

                for (var row = rowMin; row <= rowMax; row++)
                    for (var col = colMin; col <= colMax; col++)
                    {
                        var xtarget = maxOuterExtend.LngMin + row * width;
                        var ytarget = maxOuterExtend.LatMin + col * width;

                        if (polygon.Pip(xtarget, ytarget))
                        {
                            var valItem = new KigingGridValue
                            {
                                Row = row,
                                Col = col,
                                Value = Predict(new GeoPoint(xtarget, ytarget), variogram)
                            };
                            values.Add(valItem);

                            grid[row][col] = valItem.Value;
                        }

                    }
            }

            return new KigingGrid
            {
                Grid = grid,
                Values = values.ToArray(),
                Extend = maxOuterExtend,
                Width = width
            };
        }

        /// <summary>
        /// return rects that can be drawn to canvas
        /// </summary>
        /// <remarks>Only recommended to use this function when too much points to plot</remarks>
        public static PlotData[] Plot(KrigingTrainVariogram variogram, KigingGrid grid, DrawOptions options)
        {
            var result = new List<PlotData>();
            var canvas = options.Canvas;
            var extent = options.Extent;

            var extentW = extent.LngMax - extent.LngMin;
            var extentH = extent.LatMax - extent.LatMin;
            var valueDiff = variogram.MaxValue - variogram.MinValue;

            var wx = Math.Ceiling(grid.Width * canvas.Width / extentW);
            var wy = Math.Ceiling(grid.Width * canvas.Height / extentH);

            foreach (var gridVal in grid.Values)
            {
                var x = canvas.Width * (gridVal.Row * grid.Width + grid.Extend.LngMin - extent.LngMin) / extentW;
                var y = canvas.Height * (1 - (gridVal.Col * grid.Width + grid.Extend.LatMin - extent.LatMin) / extentH);

                var z = (gridVal.Value - variogram.MinValue) / valueDiff;
                z = Math.Min(1, Math.Max(0, z));

                result.Add(new PlotData
                {
                    X = Math.Round(x - wx / 2),
                    Y = Math.Round(y - wy / 2),
                    Width = wx,
                    Height = wy,
                    Color = options.Colors[(int)Math.Floor((options.Colors.Length - 1) * z)]
                });
            }

            return result.ToArray();
        }

        /// <summary>
        /// Model prediction
        /// </summary>
        public static double Predict(GeoPoint point, KrigingTrainVariogram variogram)
        {
            var temp = variogram.Data.Select(g => variogram.Model(point.Distance(g))).ToArray();
            return Matrix.Multiply(temp, variogram.M, 1, variogram.Data.Length, 1)[0];
        }

        public static double Variance(GeoPoint point, KrigingTrainVariogram variogram)
        {
            var length = variogram.Data.Length;
            var temp = variogram.Data.Select(g => variogram.Model(point.Distance(g))).ToArray();

            return variogram.Model(0) +
                Matrix.Multiply(
                    Matrix.Multiply(temp, variogram.K, 1, length, length),
                    temp, 1, length, 1)[0];
        }


        private class Matrix
        {
            public static double[] Transpose(double[] x, int n, int m)
            {
                var temp = new double[m * n];
                for (var i = 0; i < n; i++)
                    for (var j = 0; j < m; j++)
                        temp[j * n + i] = x[i * m + j];

                return temp;
            }

            /// <summary>
            /// Naive matrix multiplication
            /// </summary>
            public static double[] Multiply(double[] x, double[] y, int n, int m, int p)
            {
                var temp = new double[n * p];
                for (var i = 0; i < n; i++)
                    for (var j = 0; j < p; j++)
                    {
                        temp[i * p + j] = 0;

                        for (var k = 0; k < m; k++)
                            temp[i * p + j] += x[i * m + k] * y[k * p + j];
                    }

                return temp;
            }

            public static double[] Add(double[] x, double[] y, int n) => Add(x, y, n, n);
            public static double[] Add(double[] x, double[] y, int n, int m)
            {
                var temp = new double[n * m];
                for (var i = 0; i < n; i++)
                    for (var j = 0; j < m; j++)
                    {
                        var idx = i * m + j;

                        temp[idx] = x[idx] + y[idx];
                    }


                return temp;
            }

            /// <summary>
            /// Matrix algebra
            /// </summary>
            public static double[] Diag(double c, int n)
            {
                var temp = new double[n * n];
                for (var i = 0; i < n; i++)
                    temp[i * n + i] = c;

                return temp;
            }

            /// <summary>
            /// Cholesky decomposition
            /// </summary>
            public static bool Chol(double[] x, int n)
            {
                var temp = new double[n];
                for (var i = 0; i < n; i++)
                    temp[i] = x[i * n + i];

                for (var i = 0; i < n; i++)
                {
                    for (var j = 0; j < i; j++)
                        temp[i] -= x[i * n + j] * x[i * n + j];

                    if (temp[i] <= 0)
                        return false;

                    temp[i] = Math.Sqrt(temp[i]);
                    for (var j = i + 1; j < n; j++)
                    {
                        for (var k = 0; k < i; k++)
                            x[j * n + i] -= x[j * n + k] * x[i * n + k];

                        x[j * n + i] /= temp[i];
                    }
                }

                for (var i = 0; i < n; i++)
                    x[i * n + i] = temp[i];

                return true;
            }

            /// <summary>
            /// Inversion of cholesky decomposition
            /// </summary>
            public static void Chol2Inv(double[] x, int n)
            {
                for (var i = 0; i < n; i++)
                {
                    x[i * n + i] = 1 / x[i * n + i];

                    for (var j = i + 1; j < n; j++)
                    {
                        var sum = 0.0;
                        for (var k = i; k < j; k++)
                            sum -= x[j * n + k] * x[k * n + i];

                        x[j * n + i] = sum / x[j * n + j];
                    }
                }

                for (var i = 0; i < n; i++)
                    for (var j = i + 1; j < n; j++)
                        x[i * n + j] = 0;

                for (var i = 0; i < n; i++)
                {
                    x[i * n + i] *= x[i * n + i];

                    for (var k = i + 1; k < n; k++)
                        x[i * n + i] += x[k * n + i] * x[k * n + i];

                    for (var j = i + 1; j < n; j++)
                        for (var k = j; k < n; k++)
                            x[i * n + j] += x[k * n + i] * x[k * n + j];
                }

                for (var i = 0; i < n; i++)
                    for (var j = 0; j < i; j++)
                        x[i * n + j] = x[j * n + i];
            }

            /// <summary>
            /// Inversion via gauss-jordan elimination
            /// </summary>
            public static bool Solve(double[] x, int n)
            {
                var indxc = new int[n];
                var indxr = new int[n];
                var ipiv = new double[n];

                for (var i = 0; i < n; i++)
                {
                    var b = new double[n * n];
                    var big = 0.0;
                    int icol = 0, irow = 0;

                    for (var j = 0; j < n; j++)
                    {
                        b[i * n + j] = i == j ? 1 : 0;

                        if (ipiv[j] != 1)
                        {
                            for (var k = 0; k < n; k++)
                            {
                                if (ipiv[k] == 0)
                                {
                                    var abs = Math.Abs(x[j * n + k]);
                                    if (abs >= big)
                                    {
                                        big = abs;

                                        irow = j;
                                        icol = k;
                                    }
                                }
                            }
                        }
                    }
                    ++(ipiv[icol]);

                    if (irow != icol)
                    {
                        for (var l = 0; l < n; l++)
                        {
                            var rowIdx = irow * n + l;
                            var colIdx = icol * n + l;

                            var temp = x[rowIdx];
                            x[rowIdx] = x[colIdx];
                            x[colIdx] = temp;

                            temp = b[rowIdx];
                            b[rowIdx] = b[colIdx];
                            b[colIdx] = temp;
                        }
                    }

                    indxr[i] = irow;
                    indxc[i] = icol;

                    if (x[icol * n + icol] == 0)
                        return false; // Singular

                    var pivinv = 1 / x[icol * n + icol];
                    x[icol * n + icol] = 1;

                    for (var l = 0; l < n; l++)
                    {
                        var idx = icol * n + l;

                        x[idx] *= pivinv;
                        b[idx] *= pivinv;
                    }

                    for (var ll = 0; ll < n; ll++)
                    {
                        if (ll != icol)
                        {
                            var dum = x[ll * n + icol];
                            x[ll * n + icol] = 0;

                            for (var l = 0; l < n; l++)
                            {
                                x[ll * n + l] -= x[icol * n + l] * dum;
                                b[ll * n + l] -= b[icol * n + l] * dum;
                            }
                        }
                    }
                }

                for (var i = (n - 1); i >= 0; i--)
                    if (indxr[i] != indxc[i])
                    {
                        for (var k = 0; k < n; k++)
                        {
                            var rIdx = k * n + indxr[i];
                            var cIdx = k * n + indxc[i];

                            var temp = x[rIdx];
                            x[rIdx] = x[cIdx];
                            x[cIdx] = temp;
                        }
                    }

                return true;
            }

            public static double[] Clone(double[] from)
            {
                var temp = new double[from.Length];
                from.CopyTo(temp, 0);

                return temp;
            }
        }
    }

    public class KrigingTrainVariogram
    {
        private readonly KrigingTrainModel _model;

        public KrigingTrainVariogram(GeoDoubleData[] data, KrigingTrainModel model)
        {
            _model = model;
            Data = data;
        }

        public GeoDoubleData[] Data { get; }

        public double Nugget { get; set; }
        public double Range { get; set; }
        public double Sill { get; set; }
        public double[] K { get; set; }
        public double[] M { get; set; }

        public double MinValue { get; set; }
        public double MaxValue { get; set; }

        public double Transform(double h)
        {
            var a = h / Range;

            switch (_model)
            {
                case KrigingTrainModel.Gaussian:
                    return 1.0 - Math.Exp(-3 * Math.Pow(a, 2));
                case KrigingTrainModel.Exponential:
                    return 1.0 - Math.Exp(-3 * a);
                case KrigingTrainModel.Spherical:
                    return 1.5 * (a) - 0.5 * Math.Pow(a, 3);
            }
            return 0;
        }

        public double Model(double h)
        {
            switch (_model)
            {
                case KrigingTrainModel.Gaussian:
                    return Nugget + ((Sill - Nugget) / Range) * Transform(h);
                case KrigingTrainModel.Exponential:
                    return Nugget + ((Sill - Nugget) / Range) * Transform(h);
                case KrigingTrainModel.Spherical:
                    if (h > Range)
                        return Nugget + (Sill - Nugget) / Range;
                    return Nugget + ((Sill - Nugget) / Range) * Transform(h);
            }
            return 0;
        }
    }

    public class KigingGrid
    {
        public double?[][] Grid { get; set; }

        public KigingGridValue[] Values { get; set; }

        public GeoExtend Extend { get; set; }

        public double Width { get; set; }

        /// <summary>
        /// convert to  return type of [kriging.js] gird function
        /// </summary>
        /// <returns></returns>
        public object ForPlot(KrigingTrainVariogram variogram)
        {
            dynamic result = new ExpandoObject();

            IDictionary<string, object> values = result;
            for (int i = 0; i < Grid.Length; i++)
                values[i.ToString()] = Grid[i];

            result.length = Grid.Length;
            result.width = Width;
            result.xlim = new double[2] { Extend.LngMin, Extend.LngMax };
            result.ylim = new double[2] { Extend.LatMin, Extend.LatMax };
            result.zlim = new double[2] { variogram.MinValue, variogram.MaxValue };

            return result;
        }
    }
    [DebuggerDisplay("[{Row},{Col}]{Value}")]
    public class KigingGridValue
    {
        public int Row { get; set; }
        public int Col { get; set; }
        public double Value { get; set; }
    }


    public class DrawOptions
    {
        public SizeInfo Canvas { get; set; }
        public GeoExtend Extent { get; set; }

        public string[] Colors { get; set; }
    }
    public class SizeInfo
    {
        public double Width { get; set; }
        public double Height { get; set; }
    }
    [DebuggerDisplay("({X},{Y}){Color}")]
    public class PlotData : SizeInfo
    {
        public double X { get; set; }
        public double Y { get; set; }

        public string Color { get; set; }
    }
}
