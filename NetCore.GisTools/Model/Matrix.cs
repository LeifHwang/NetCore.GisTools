using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;

namespace NetCore.GisTools.Model
{
    [DebuggerDisplay("Matrix ({RowCount} * {ColumnCount})")]
    [DebuggerTypeProxy(typeof(MatrixDebugView))]
    public sealed class Matrix : IEnumerable<double>
    {
        private readonly double[,] _data;

        public Matrix(int rows, int columns)
        {
            _data = new double[rows, columns];

            RowCount = rows;
            ColumnCount = columns;
        }
        /// <summary>
        /// square matrix
        /// </summary>
        /// <param name="rows"></param>
        public Matrix(int rows): this(rows,rows) { }

        /// <summary>
        /// clone
        /// </summary>
        /// <param name="matrix"></param>
        public Matrix(Matrix matrix)
        {
            RowCount = matrix._data.GetUpperBound(0) + 1;
            ColumnCount = matrix._data.GetUpperBound(1) + 1;

            _data = new double[RowCount, ColumnCount];
            Array.Copy(matrix._data, _data, RowCount * ColumnCount);
        }
        public Matrix(IEnumerable<double> values)
        {
            RowCount = values.Count();
            ColumnCount = 1;

            var index = 0;
            _data = new double[RowCount, ColumnCount];
            foreach (var item in values)
            {
                _data[index++, 0] = item;
            }
        }

        public int RowCount { get; }
        public int ColumnCount { get; }
        public double this[int row, int column]
        {
            get { return _data[row, column]; }
            set { _data[row, column] = value; }
        }

        public IEnumerator<double> GetEnumerator()
        {
            foreach (var item in _data)
                yield return item;
        }
        IEnumerator IEnumerable.GetEnumerator()
        {
            yield return this.AsEnumerable();
        }


        /// <summary>
        /// create a matrix with values on the principal diagonal
        /// </summary>
        /// <param name="value"></param>
        /// <param name="size"></param>
        /// <returns></returns>
        public static Matrix Diagonal(double value, int size)
        {
            var result = new Matrix(size);
            for (var i = 0; i < size; i++)
                result[i, i] = value;
            return result;
        }

        public static Matrix Transpose(Matrix matrix)
        {
            var result = new Matrix(matrix.ColumnCount, matrix.RowCount);
            for (int r = 0; r < matrix.RowCount; r++)
            {
                for (int c = 0; c < matrix.ColumnCount; c++)
                    result[c, r] = matrix[r, c];
            }
            return result;
        }

        public static Matrix operator +(Matrix a, Matrix b)
        {
            if (a.RowCount != b.RowCount || a.ColumnCount != b.ColumnCount)
                throw new Exception("Matrix dimension mismatching!");

            var result = new Matrix(a.RowCount, a.ColumnCount);
            for (int i = 0; i < a.RowCount; i++)
                for (int j = 0; j < a.ColumnCount; j++)
                    result[i, j] = a[i, j] + b[i, j];
            return result;
        }
        public static Matrix operator *(Matrix a, Matrix b)
        {
            //if (a.ColumnCount != b.RowCount)
            //    throw new Exception("Matrix dimension mismatching!");

            var result = new Matrix(a.RowCount, b.ColumnCount);
            for (int i = 0; i < a.RowCount; i++)
            {
                for (int j = 0; j < b.ColumnCount; j++)
                {
                    for (int k = 0; k < a.ColumnCount; k++)
                        result[i, j] += a[i, k] * b[k, j];
                }
            }
            return result;
        }
    }

    class MatrixDebugView
    {
        [DebuggerBrowsable(DebuggerBrowsableState.Never)]
        private readonly Matrix _matrix;

        public MatrixDebugView(Matrix matrix) => _matrix = matrix;

        [DebuggerBrowsable(DebuggerBrowsableState.RootHidden)]
        private string[] Rows
        {
            get
            {
                var temp = new List<string>();
                for (int i = 0; i < _matrix.RowCount; i++)
                {
                    var row = new string[_matrix.ColumnCount];
                    for (int j = 0; j < _matrix.ColumnCount; j++)
                    {
                        row[j] = _matrix[i, j].ToString();
                    }
                    temp.Add(string.Join(", ", row));
                }
                return temp.ToArray();
            }
        }
    }
}
