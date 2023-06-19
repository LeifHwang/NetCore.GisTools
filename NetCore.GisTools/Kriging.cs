using NetCore.GisTools.Model;

using System;
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

            var result = new KrigingTrainVariogram(values, model);

            // Lag distance/semivariance
            var inputPointCount = values.Length;
            var arithmeticSum = (inputPointCount * inputPointCount - inputPointCount) / 2;
            var distance = new double[arithmeticSum][];
            var inputValues = new double[inputPointCount];
            for (int i = 0, k = 0; i < inputPointCount; i++)
            {
                inputValues[i] = values[i].Value;
                for (int j = 0; j < i; j++, k++)
                {
                    var temp = new double[2];
                    temp[0] = values[i].Distance(values[j]);
                    temp[1] = Math.Abs(values[i].Value - values[j].Value);

                    distance[k] = temp;
                }
            }
            Array.Sort(distance, (a, b) => a[0].CompareTo(b[0]));

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

            var cloneLaglagModelRanges = new double[lagModelRanges.Length];
            lagModelRanges.CopyTo(cloneLaglagModelRanges, 0);

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
            var cloneSigmaModels = new double[sigmaModels.Length];
            sigmaModels.CopyTo(cloneSigmaModels, 0);

            if (Matrix.Chol(sigmaModels, inputPointCount))
                Matrix.Chol2Inv(sigmaModels, inputPointCount);
            else
            {
                Matrix.Solve(cloneSigmaModels, inputPointCount);
                sigmaModels = cloneSigmaModels;
            }

            // Copy unprojected inverted matrix as K
            result.K = new double[sigmaModels.Length];
            sigmaModels.CopyTo(result.K, 0);

            result.M = Matrix.Multiply(sigmaModels, inputValues, inputPointCount, inputPointCount, 1);

            return result;
        }

        /// <summary>
        /// Gridded matrices or contour paths
        /// </summary>
        public static KigingGrid Grid(GeoPolygon[] polygons, KrigingTrainVariogram variogram, double width)
        {
            if (polygons.Length == 0)
                return new KigingGrid();

            var filteredPolygons = polygons.Where(p => p.Points.Length > 0).ToArray();
            if (filteredPolygons.Length == 0)
                return new KigingGrid();

            // Boundaries of polygons space
            var xlim = new double[2] { filteredPolygons[0].Points[0].Longitude, filteredPolygons[0].Points[0].Longitude };
            var ylim = new double[2] { filteredPolygons[0].Points[0].Latitude, filteredPolygons[0].Points[0].Latitude };
            for (var i = 0; i < filteredPolygons.Length; i++)
                for (var j = 0; j < filteredPolygons[i].Points.Length; j++)
                {
                    var longitude = filteredPolygons[i].Points[j].Longitude;
                    var latitude = filteredPolygons[i].Points[j].Latitude;

                    if (longitude < xlim[0])
                        xlim[0] = longitude;
                    if (longitude > xlim[1])
                        xlim[1] = longitude;

                    if (latitude < ylim[0])
                        ylim[0] = latitude;
                    if (latitude > ylim[1])
                        ylim[1] = latitude;
                }

            // Alloc for O(n^2) space
            var lxlim = new double[2];
            var lylim = new double[2];

            var x = (int)Math.Ceiling((xlim[1] - xlim[0]) / width);
            var y = (int)Math.Ceiling((ylim[1] - ylim[0]) / width);

            var grid = new double[x + 1][];
            for (var i = 0; i <= x; i++)
                grid[i] = new double[y + 1];

            for (var i = 0; i < polygons.Length; i++)
            {
                // Range for polygons[i]
                lxlim[0] = filteredPolygons[i].Points[0].Longitude;
                lxlim[1] = lxlim[0];
                lylim[0] = filteredPolygons[i].Points[0].Latitude;
                lylim[1] = lylim[0];

                for (var j = 1; j < filteredPolygons[i].Points.Length; j++)
                {
                    var longitude = filteredPolygons[i].Points[j].Longitude;
                    var latitude = filteredPolygons[i].Points[j].Latitude;

                    if (longitude < lxlim[0])
                        lxlim[0] = longitude;
                    if (longitude > lxlim[1])
                        lxlim[1] = longitude;

                    if (latitude < lylim[0])
                        lylim[0] = latitude;
                    if (latitude > lylim[1])
                        lylim[1] = latitude;
                }

                // Loop through polygon subspace
                var a = new int[2];
                a[0] = (int)Math.Floor(((lxlim[0] - ((lxlim[0] - xlim[0]) % width)) - xlim[0]) / width);
                a[1] = (int)Math.Ceiling(((lxlim[1] - ((lxlim[1] - xlim[1]) % width)) - xlim[0]) / width);

                var b = new int[2];
                b[0] = (int)Math.Floor(((lylim[0] - ((lylim[0] - ylim[0]) % width)) - ylim[0]) / width);
                b[1] = (int)Math.Ceiling(((lylim[1] - ((lylim[1] - ylim[1]) % width)) - ylim[0]) / width);

                for (var j = a[0]; j <= a[1]; j++)
                    for (var k = b[0]; k <= b[1]; k++)
                    {
                        var xtarget = xlim[0] + j * width;
                        var ytarget = ylim[0] + k * width;
                        if (filteredPolygons[i].Pip(xtarget, ytarget))
                            grid[j][k] = Predict(new LocationPoint(xtarget, ytarget), variogram);
                    }
            }

            return new KigingGrid
            {
                Grid = grid,
                Xlim = xlim,
                Ylim = ylim,
                Width = width
            };
        }

        /// <summary>
        /// Model prediction
        /// </summary>
        public static double Predict(LocationPoint point, KrigingTrainVariogram variogram)
        {
            var temp = variogram.Data.Select(g => variogram.Model(point.Distance(g))).ToArray();
            return Matrix.Multiply(temp, variogram.M, 1, variogram.Data.Length, 1)[0];
        }

        public static double Variance(LocationPoint point, KrigingTrainVariogram variogram)
        {
            var length = variogram.Data.Length;
            var temp = variogram.Data.Select(g => variogram.Model(point.Distance(g))).ToArray();

            return variogram.Model(0) +
                Matrix.Multiply(
                    Matrix.Multiply(temp, variogram.K, 1, length, length),
                    temp, 1, length, 1)[0];
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
        public double[][] Grid { get; set; }

        public double Width { get; set; }

        public double[] Xlim { get; set; }
        public double[] Ylim { get; set; }

        //public double[] Zlim { get; set; }
    }

    internal class Matrix
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
            {
                for (var j = 0; j < p; j++)
                {
                    temp[i * p + j] = 0;

                    for (var k = 0; k < m; k++)
                        temp[i * p + j] += x[i * m + k] * y[k * p + j];
                }
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

            for (var l = (n - 1); l >= 0; l--)
                if (indxr[l] != indxc[l])
                {
                    for (var k = 0; k < n; k++)
                    {
                        var rIdx = k * n + indxr[l];
                        var cIdx = k * n + indxc[l];

                        var temp = x[rIdx];
                        x[rIdx] = x[cIdx];
                        x[cIdx] = temp;
                    }
                }

            return true;
        }
    }
}
