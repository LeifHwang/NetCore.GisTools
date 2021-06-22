using NetCore.GisTools.Model;
using NetCore.GisTools.Utils;

using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace NetCore.GisTools
{
    public static class Kriging
    {
        public async static Task<KrigingTrainVariogram> TrainAsync(KrigingInputModel input, int sigma, int alpha,
            KrigingTrainModel model = KrigingTrainModel.Exponential)
        {
            var inputPointCount = input.T.Length;
            var result = new KrigingTrainVariogram(model) { A = 1m / 3, N = inputPointCount };

            // Lag distance/semi-variance
            var distance = await Task.Run(() =>
            {
                var temp = new List<double[]>(((inputPointCount * inputPointCount) - inputPointCount) / 2);
                for (int i = 0; i < inputPointCount; i++)
                {
                    for (int j = 0; j < i; j++)
                    {
                        var dis = GeoMath.Distance(input.X[i], input.Y[i], input.X[j], input.Y[j]);
                        temp.Add(new double[] { dis, Math.Abs(input.T[i] - input.T[j]) });
                    }
                }
                temp.Sort((a, b) => a[0].CompareTo(b[0]));
                return temp;
            });
            result.Range = distance[^1][0];

            // Bin lag distance
            var lagsCount = Math.Min(distance.Count, 30);
            var tolerance = result.Range / lagsCount;
            var lags = new double[lagsCount];
            var semis = new double[lagsCount];
            if (lagsCount < 30)
            {
                for (var l = 0; l < lagsCount; l++)
                {
                    lags[l] = distance[l][0];
                    semis[l] = distance[l][1];
                }
            }
            else
            {
                int l = 0;
                for (int i = 0, j = 0, k = 0; i < lagsCount && j < distance.Count; i++, k = 0)
                {
                    while (distance[j][0] <= ((i + 1) * tolerance))
                    {
                        lags[l] += distance[j][0];
                        semis[l] += distance[j][1];

                        j++; k++;
                        if (j >= distance.Count) break;
                    }

                    if (k > 0)
                    {
                        lags[l] /= k;
                        semis[l] /= k;
                        l++;
                    }
                }
                if (l < 2) return result; // Not enough points
            }

            result.Range = lags[lagsCount - 1] - lags[0];

            // Feature transformation
            var lagModels = new double[lagsCount * 2];
            for (var i = 0; i < lagsCount; i++)
            {
                lagModels[i * 2] = 1.0;
                lagModels[(i * 2) + 1] = result.Model1(lags[i]);
            }

            // Least squares
            var lagModelsTrans = MatrixTranspose(lagModels, lagsCount, 2);
            var lagModelRanges = MatrixMultiply(lagModelsTrans, lagModels, 2, lagsCount, 2);
            lagModelRanges = MatrixAdd(lagModelRanges, MatrixDiag(1.0 / alpha, 2), 2, 2);

            var cloneLaglagModelRanges = new double[lagModelRanges.Length];
            lagModelRanges.CopyTo(cloneLaglagModelRanges, 0);

            if (MatrixCholesky(lagModelRanges, 2))
            {
                MatrixCholesky2Inversion(lagModelRanges, 2);
            }
            else
            {
                MatrixSolve(lagModels, 2);
                lagModelRanges = cloneLaglagModelRanges;
            }

            var rectangle = MatrixMultiply(MatrixMultiply(lagModelRanges, lagModelsTrans, 2, 2, lagsCount), semis, 2, lagsCount, 1);

            // Variogram parameters
            result.Nugget = rectangle[0];
            result.Sill = (rectangle[1] * result.Range) + result.Nugget;
            result.N = inputPointCount;

            // Gram matrix with prior
            var models = new double[inputPointCount * inputPointCount];
            Parallel.For(0, inputPointCount, i =>
            {
                for (var j = 0; j < i; j++)
                {
                    var model = result.Model(GeoMath.Distance(input.X[i], input.Y[i], input.X[j], input.Y[j]));

                    models[(i * inputPointCount) + j] = model;
                    models[(j * inputPointCount) + i] = model;
                }
                models[(i * inputPointCount) + i] = result.Model(0);
            });

            // Inverse penalized Gram matrix projected to target vector
            var sigmaModels = MatrixAdd(models, MatrixDiag(sigma, inputPointCount), inputPointCount, inputPointCount);
            var cloneSigmaModels = new double[sigmaModels.Length];
            sigmaModels.CopyTo(cloneSigmaModels, 0);

            if (MatrixCholesky(sigmaModels, inputPointCount))
            {
                MatrixCholesky2Inversion(sigmaModels, inputPointCount);
            }
            else
            {
                MatrixSolve(lagModels, 2);
                sigmaModels = cloneSigmaModels;
            }

            // Copy unprojected inverted matrix as K
            result.K = new double[sigmaModels.Length];
            sigmaModels.CopyTo(result.K, 0);

            result.M = MatrixMultiply(sigmaModels, input.T, inputPointCount, inputPointCount, 1);

            return result;
        }

        public static KigingGrid Grid(KrigingInputModel model, KrigingTrainVariogram variogram, double width)
        {
            var result = new KigingGrid();

            var sortValues = model.T.OrderBy(a => a).ToArray();
            result.Zlim = new double[2] { sortValues[0], sortValues[sortValues.Length - 1] };

            // Boundaries of polygons space
            var xLim = new double[2] { model.MinX, model.MaxX };
            var yLim = new double[2] { model.MinY, model.MaxY };

            var polygon = new Polygon(new double[5][]
            {
                new double[2] { model.MinX, model.MinY },
                new double[2] { model.MinX, model.MaxY },
                new double[2] { model.MaxX, model.MaxY },
                new double[2] { model.MaxX, model.MinY },
                new double[2] { model.MinX, model.MinY }
            });

            var x = (int)Math.Ceiling((xLim[1] - xLim[0]) / width);
            var y = (int)Math.Ceiling((yLim[1] - yLim[0]) / width);

            // Loop through polygon subspace
            var grid = new double?[x + 1][];
            Parallel.For(0, x + 1, i =>
            {
                var index = i;
                var xTarget = xLim[0] + (index * width);

                grid[index] = new double?[y + 1];
                for (var j = 0; j <= y; j++)
                {
                    var yTarget = yLim[0] + (j * width);

                    if (polygon.Pip(xTarget, yTarget))
                    {
                        var val = 0.0;
                        for (var k = 0; k < variogram.N; k++)
                        {
                            var dis = GeoMath.Distance(model.X[k], model.Y[k], xTarget, yTarget);
                            val += variogram.Model(dis) * variogram.M[k];
                        }
                        grid[index][j] = val;
                    }
                }
            });

            result.Grid = grid;
            result.Xlim = xLim;
            result.Ylim = yLim;
            result.Width = width;

            return result;
        }


        private static double[] MatrixTranspose(double[] x, int n, int m)
        {
            var temp = new double[n * m];
            Parallel.For(0, n, i =>
            {
                for (var j = 0; j < m; j++)
                    temp[(j * n) + i] = x[(i * m) + j];
            });
            return temp;
        }
        private static double[] MatrixMultiply(double[] x, double[] y, int n, int m, int p)
        {
            var temp = new double[n * p];
            Parallel.For(0, n, (i) =>
            {
                for (var j = 0; j < p; j++)
                {
                    temp[(i * p) + j] = 0;

                    for (var k = 0; k < m; k++)
                        temp[(i * p) + j] += x[(i * m) + k] * y[(k * p) + j];
                }
            });
            return temp;
        }
        private static double[] MatrixAdd(double[] x, double[] y, int n, int m)
        {
            var temp = new double[n * m];
            Parallel.For(0, n, i =>
            {
                for (var j = 0; j < m; j++)
                    temp[(i * m) + j] = x[(i * m) + j] + y[(i * m) + j];
            });
            return temp;
        }
        // Matrix algebra
        private static double[] MatrixDiag(double c, int n)
        {
            var temp = new double[n * n];
            Parallel.For(0, n, (i) => temp[(i * n) + i] = c);
            return temp;
        }
        // Cholesky decomposition
        private static bool MatrixCholesky(double[] x, int n)
        {
            //var temp = new double[n];

            //Parallel.For(0, n, i => temp[i] = x[(i * n) + i]);

            for (var i = 0; i < n; i++)
            {
                var index = (i * n) + i;
                var temp = x[index];

                var result = Parallel.For(0, i, (j, s) =>
                {
                    var item = x[(i * n) + j];
                    temp -= item * item;

                    if (temp <= 0)
                        s.Break();
                });

                if (!result.IsCompleted)
                    return false;

                temp = Math.Sqrt(temp);

                Parallel.For(i + 1, n, j =>
                {
                    var sum = 0.0;
                    for (var k = 0; k < i; k++)
                        sum += x[(j * n) + k] * x[(i * n) + k];

                    var idx = (j * n) + i;
                    x[idx] = (x[idx] - sum) / temp;
                });

                x[index] = temp;
            }

            //Parallel.For(0, n, (i) => x[(i * n) + i] = temp[i]);
            return true;
        }
        // Inversion of cholesky decomposition
        private static void MatrixCholesky2Inversion(double[] x, int n)
        {
            for (var i = 0; i < n; i++)
            {
                var index = (i * n) + i;
                x[index] = 1.0 / x[index];

                for (var j = i + 1; j < n; j++)
                {
                    double sum = 0;
                    for (var k = i; k < j; k++)
                        sum -= x[(j * n) + k] * x[(k * n) + i];

                    x[(j * n) + i] = sum / x[(j * n) + j];
                    x[(i * n) + j] = 0;
                }
            }

            for (var i = 0; i < n; i++)
            {
                var index = (i * n) + i;
                x[index] *= x[index];

                for (var j = i + 1; j < n; j++)
                {
                    var item = x[(j * n) + i];
                    x[index] += item * item;

                    for (var k = j; k < n; k++)
                        x[(i * n) + j] += x[(k * n) + i] * x[(k * n) + j];
                }
            }

            for (var i = 0; i < n; i++)
            {
                for (var j = 0; j < i; j++)
                    x[(i * n) + j] = x[(j * n) + i];
            }
        }
        // Inversion via gauss-jordan elimination
        private static bool MatrixSolve(double[] x, int n)
        {
            var m = n;
            var b = new double[n * n];
            var ipiv = new double[n];
            var indxc = new int[n];
            var indxr = new int[n];
            int icol = 0, irow = 0;

            for (var i = 0; i < n; i++)
                for (var j = 0; j < n; j++)
                {
                    if (i == j)
                        b[(i * n) + j] = 1;
                    else
                        b[(i * n) + j] = 0;
                }

            for (var j = 0; j < n; j++)
                ipiv[j] = 0;

            for (var i = 0; i < n; i++)
            {
                double big = 0;
                for (var j = 0; j < n; j++)
                {
                    if (ipiv[j] != 1)
                    {
                        for (var k = 0; k < n; k++)
                        {
                            if (ipiv[k] == 0)
                            {
                                if (Math.Abs(x[(j * n) + k]) >= big)
                                {
                                    big = Math.Abs(x[(j * n) + k]);
                                    irow = j;
                                    icol = k;
                                }
                            }
                        }
                    }
                }
                ++ipiv[icol];

                if (irow != icol)
                {
                    for (var l = 0; l < n; l++)
                    {
                        var temp = x[(irow * n) + l];
                        x[(irow * n) + l] = x[(icol * n) + l];
                        x[(icol * n) + l] = temp;
                    }
                    for (var l = 0; l < m; l++)
                    {
                        var temp = b[(irow * n) + l];
                        b[(irow * n) + l] = b[(icol * n) + l];
                        b[(icol * n) + l] = temp;
                    }
                }
                indxr[i] = irow;
                indxc[i] = icol;

                if (x[(icol * n) + icol] == 0)
                    return false; // Singular

                var pivinv = 1 / x[(icol * n) + icol];
                x[(icol * n) + icol] = 1;

                for (var l = 0; l < n; l++)
                    x[(icol * n) + l] *= pivinv;
                for (var l = 0; l < m; l++)
                    b[(icol * n) + l] *= pivinv;

                for (var ll = 0; ll < n; ll++)
                {
                    if (ll != icol)
                    {
                        var dum = x[(ll * n) + icol];
                        x[(ll * n) + icol] = 0;

                        for (var l = 0; l < n; l++)
                            x[(ll * n) + l] -= x[(icol * n) + l] * dum;

                        for (var l = 0; l < m; l++)
                            b[(ll * n) + l] -= b[(icol * n) + l] * dum;
                    }
                }
            }
            for (var l = n - 1; l >= 0; l--)
                if (indxr[l] != indxc[l])
                {
                    for (var k = 0; k < n; k++)
                    {
                        var temp = x[(k * n) + indxr[l]];
                        x[(k * n) + indxr[l]] = x[(k * n) + indxc[l]];
                        x[(k * n) + indxc[l]] = temp;
                    }
                }
            return true;
        }
    }

    public enum KrigingTrainModel
    {
        Gaussian,
        Exponential,
        Spherical
    }

    public class KrigingTrainVariogram
    {
        private readonly KrigingTrainModel _model;

        private decimal _a;
        private double _tansA;

        private double _sill;
        private double _nugget;
        private double _range;
        private double _snDividedByR;

        public KrigingTrainVariogram(KrigingTrainModel model)
        {
            _model = model;
        }

        public double Nugget
        {
            get => _nugget; set
            {
                _nugget = value;
                _snDividedByR = (Sill - value) / Range;
            }
        }
        public double Range
        {
            get => _range; set
            {
                _range = value;
                _snDividedByR = (Sill - Nugget) / value;
            }
        }
        public double Sill
        {
            get => _sill; set
            {
                _sill = value;
                _snDividedByR = (value - Nugget) / Range;
            }
        }

        public decimal A
        {
            get => _a;
            set
            {
                _a = value;
                _tansA = (double)(1m / A);
            }
        }
        public double[] K { get; set; }
        public double[] M { get; set; }
        public int N { get; set; }

        public double Model1(double h)
        {
            h /= Range;

            //var temp = Math.Exp(-_tansA * h);

            switch (_model)
            {
                case KrigingTrainModel.Gaussian:
                    return 1.0 - Math.Exp(-_tansA * h * h);
                case KrigingTrainModel.Exponential:
                    return 1.0 - Math.Exp(-_tansA * h);
                case KrigingTrainModel.Spherical:
                    return (1.5 * h) - (0.5 * h * h * h);
            }
            return 0;
        }
        public double Model(double h)
        {
            if (h == 0)
                return Nugget;

            switch (_model)
            {
                case KrigingTrainModel.Gaussian:
                    return Nugget + (_snDividedByR * Model1(h));
                case KrigingTrainModel.Exponential:
                    return Nugget + (_snDividedByR * Model1(h));
                case KrigingTrainModel.Spherical:
                    if (h > Range)
                        return Nugget + _snDividedByR;
                    return Nugget + (_snDividedByR * Model1(h));
            }
            return 0;
        }

        //static double Exponential(double x)
        //{
        //    double sum = 1;
        //    for (int i = 15 - 1; i > 0; --i)
        //        sum = 1 + (x * sum / i);
        //    return sum;
        //}
    }

    public class KigingGrid
    {
        public double?[][] Grid { get; set; }

        public double Width { get; set; }

        public double[] Xlim { get; set; }
        public double[] Ylim { get; set; }
        public double[] Zlim { get; set; }
    }

    class Polygon
    {
        private readonly double[][] _points;

        public Polygon(double[][] points)
        {
            _points = points;
        }

        public bool Pip(double x, double y)
        {
            for (int i = 0, j = _points.Length - 1; i < _points.Length; j = i++)
            {
                var current = _points[i];
                var previous = _points[j];

                if (((current[1] > y) != (previous[1] > y)) &&
                    (x < ((previous[0] - current[0]) * (y - current[1]) / (previous[1] - current[1])) + current[0]))
                {
                    return true;
                }
            }
            return false;
        }
    }
}
