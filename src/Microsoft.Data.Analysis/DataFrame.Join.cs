// Licensed to the .NET Foundation under one or more agreements.
// The .NET Foundation licenses this file to you under the MIT license.
// See the LICENSE file in the project root for more information.

using System;
using System.Collections.Generic;
using System.Linq;

namespace Microsoft.Data.Analysis
{
    public enum JoinAlgorithm
    {
        Left,
        Right,
        FullOuter,
        Inner
    }

    /// <summary>
    /// A DataFrame to support indexing, binary operations, sorting, selection and other APIs. This will eventually also expose an IDataView for ML.NET
    /// </summary>
    public partial class DataFrame
    {

        private void SetSuffixForDuplicatedColumnNames(DataFrame dataFrame, DataFrameColumn column, string leftSuffix, string rightSuffix)
        {
            int index = dataFrame._columnCollection.IndexOf(column.Name);
            while (index != -1)
            {
                // Pre-existing column. Change name
                DataFrameColumn existingColumn = dataFrame.Columns[index];
                dataFrame._columnCollection.SetColumnName(existingColumn, existingColumn.Name + leftSuffix);
                column.SetName(column.Name + rightSuffix);
                index = dataFrame._columnCollection.IndexOf(column.Name);
            }
        }

        /// <summary>
        /// Joins columns of another <see cref="DataFrame"/>
        /// </summary>
        /// <param name="other">The other <see cref="DataFrame"/> to join.</param>
        /// <param name="leftSuffix">The suffix to add to this <see cref="DataFrame"/>'s column if there are common column names</param>
        /// <param name="rightSuffix">The suffix to add to the <paramref name="other"/>'s column if there are common column names</param>
        /// <param name="joinAlgorithm">The <see cref="JoinAlgorithm"/> to use.</param>
        /// <returns>A new <see cref="DataFrame"/></returns>
        public DataFrame Join(DataFrame other, string leftSuffix = "_left", string rightSuffix = "_right", JoinAlgorithm joinAlgorithm = JoinAlgorithm.Left)
        {
            DataFrame ret = new DataFrame();
            if (joinAlgorithm == JoinAlgorithm.Left)
            {
                for (int i = 0; i < Columns.Count; i++)
                {
                    DataFrameColumn newColumn = Columns[i].Clone();
                    ret.Columns.Insert(ret.Columns.Count, newColumn);
                }
                long minLength = Math.Min(Rows.Count, other.Rows.Count);
                PrimitiveDataFrameColumn<long> mapIndices = new PrimitiveDataFrameColumn<long>("mapIndices", minLength);
                for (long i = 0; i < minLength; i++)
                {
                    mapIndices[i] = i;
                }
                for (int i = 0; i < other.Columns.Count; i++)
                {
                    DataFrameColumn newColumn;
                    if (other.Rows.Count < Rows.Count)
                    {
                        newColumn = other.Columns[i].Clone(numberOfNullsToAppend: Rows.Count - other.Rows.Count);
                    }
                    else
                    {
                        newColumn = other.Columns[i].Clone(mapIndices);
                    }
                    SetSuffixForDuplicatedColumnNames(ret, newColumn, leftSuffix, rightSuffix);
                    ret.Columns.Insert(ret.Columns.Count, newColumn);
                }
            }
            else if (joinAlgorithm == JoinAlgorithm.Right)
            {
                long minLength = Math.Min(Rows.Count, other.Rows.Count);
                PrimitiveDataFrameColumn<long> mapIndices = new PrimitiveDataFrameColumn<long>("mapIndices", minLength);
                for (long i = 0; i < minLength; i++)
                {
                    mapIndices[i] = i;
                }
                for (int i = 0; i < Columns.Count; i++)
                {
                    DataFrameColumn newColumn;
                    if (Rows.Count < other.Rows.Count)
                    {
                        newColumn = Columns[i].Clone(numberOfNullsToAppend: other.Rows.Count - Rows.Count);
                    }
                    else
                    {
                        newColumn = Columns[i].Clone(mapIndices);
                    }
                    ret.Columns.Insert(ret.Columns.Count, newColumn);
                }
                for (int i = 0; i < other.Columns.Count; i++)
                {
                    DataFrameColumn newColumn = other.Columns[i].Clone();
                    SetSuffixForDuplicatedColumnNames(ret, newColumn, leftSuffix, rightSuffix);
                    ret.Columns.Insert(ret.Columns.Count, newColumn);
                }
            }
            else if (joinAlgorithm == JoinAlgorithm.FullOuter)
            {
                long newRowCount = Math.Max(Rows.Count, other.Rows.Count);
                long numberOfNulls = newRowCount - Rows.Count;
                for (int i = 0; i < Columns.Count; i++)
                {
                    DataFrameColumn newColumn = Columns[i].Clone(numberOfNullsToAppend: numberOfNulls);
                    ret.Columns.Insert(ret.Columns.Count, newColumn);
                }
                numberOfNulls = newRowCount - other.Rows.Count;
                for (int i = 0; i < other.Columns.Count; i++)
                {
                    DataFrameColumn newColumn = other.Columns[i].Clone(numberOfNullsToAppend: numberOfNulls);
                    SetSuffixForDuplicatedColumnNames(ret, newColumn, leftSuffix, rightSuffix);
                    ret.Columns.Insert(ret.Columns.Count, newColumn);
                }
            }
            else if (joinAlgorithm == JoinAlgorithm.Inner)
            {
                long newRowCount = Math.Min(Rows.Count, other.Rows.Count);
                PrimitiveDataFrameColumn<long> mapIndices = new PrimitiveDataFrameColumn<long>("mapIndices", newRowCount);
                for (long i = 0; i < newRowCount; i++)
                {
                    mapIndices[i] = i;
                }
                for (int i = 0; i < Columns.Count; i++)
                {
                    DataFrameColumn newColumn = Columns[i].Clone(mapIndices);
                    ret.Columns.Insert(ret.Columns.Count, newColumn);
                }
                for (int i = 0; i < other.Columns.Count; i++)
                {
                    DataFrameColumn newColumn = other.Columns[i].Clone(mapIndices);
                    SetSuffixForDuplicatedColumnNames(ret, newColumn, leftSuffix, rightSuffix);
                    ret.Columns.Insert(ret.Columns.Count, newColumn);
                }
            }
            return ret;
        }

        private static bool IsAnyNullValueInColumns(IReadOnlyCollection<DataFrameColumn> columns, long index)
        {
            foreach (var column in columns)
            {
                if (column[index] == null)
                    return true;
            }
            return false;
        }

        /// <summary> 
        /// Merge DataFrames with a database style join (for backward compatibility)
        /// </summary> 
        /// <param name="other"></param> 
        /// <param name="leftJoinColumn"></param> 
        /// <param name="rightJoinColumn"></param> 
        /// <param name="leftSuffix"></param> 
        /// <param name="rightSuffix"></param> 
        /// <param name="joinAlgorithm"></param> 
        /// <param name="leaveColumns"></param> 
        /// <param name="shouldAcceptRow"></param> 
        /// <returns></returns> 
        public DataFrame Merge<TKey>(DataFrame other, string leftJoinColumn, string rightJoinColumn, string leftSuffix = "_left", string rightSuffix = "_right",
            JoinAlgorithm joinAlgorithm = JoinAlgorithm.Left, IReadOnlyCollection<string> leaveColumns = null, Func<long?, long?, bool> shouldAcceptRow = null)
        {
            return Merge(other, new[] { leftJoinColumn }, new[] { rightJoinColumn }, leftSuffix, rightSuffix,
                joinAlgorithm, leaveColumns, shouldAcceptRow);
        }

        private static HashSet<long> Merge(DataFrame retainedDataFrame, DataFrame supplementaryDataFrame,
            string[] retainedJoinColumnNames, string[] supplemetaryJoinColumnNames,
            out PrimitiveDataFrameColumn<long> retainedRowIndices, out PrimitiveDataFrameColumn<long> supplementaryRowIndices,
            bool isLeftDataFrameRetained, bool isInner = false, bool calculateIntersection = false, Func<long?, long?, bool> shouldAcceptRow = null)
        {
            if (retainedJoinColumnNames == null)
                throw new ArgumentNullException(nameof(retainedJoinColumnNames));

            if (supplemetaryJoinColumnNames == null)
                throw new ArgumentNullException(nameof(supplemetaryJoinColumnNames));

            if (retainedJoinColumnNames.Length != supplemetaryJoinColumnNames.Length)
                throw new ArgumentException(Strings.MismatchedArrayLengths, nameof(retainedJoinColumnNames));

            Dictionary<long, ICollection<long>> occurrences = GetOccurences(retainedDataFrame, supplementaryDataFrame,
                retainedJoinColumnNames, supplemetaryJoinColumnNames, out HashSet<long> supplementaryJoinColumnsNullIndices);

            return PerformMerging(retainedDataFrame, retainedJoinColumnNames, occurrences, supplementaryJoinColumnsNullIndices,
                out retainedRowIndices, out supplementaryRowIndices, isLeftDataFrameRetained, isInner, calculateIntersection, shouldAcceptRow);
        }

        private static Dictionary<long, ICollection<long>> GetOccurences(DataFrame retainedDataFrame, DataFrame supplementaryDataFrame,
            string[] retainedJoinColumnNames, string[] supplemetaryJoinColumnNames, out HashSet<long> supplementaryJoinColumnsNullIndices)
        {
            supplementaryJoinColumnsNullIndices = new HashSet<long>();

            // Get occurrences of values in columns used for join in the retained and supplementary dataframes

            Dictionary<long, ICollection<long>> occurrences = null;
            Dictionary<long, long> retainedIndicesReverseMapping = null;

            for (int colNameIndex = 0; colNameIndex < retainedJoinColumnNames.Length; colNameIndex++)
            {
                DataFrameColumn shrinkedRetainedColumn = retainedDataFrame.Columns[retainedJoinColumnNames[colNameIndex]];

                // Shrink retained column by row occurrences from previous step
                if (occurrences != null)
                {
                    // Only rows with occurences from previose step should go for futher processing
                    var shrinkedRetainedIndices = occurrences.Keys.ToArray();

                    // Create reverse mapping of index of the row in the shrinked column to the index of this row in the original dataframe (new index -> original index)
                    var newRetainedIndicesReverseMapping = new Dictionary<long, long>(shrinkedRetainedIndices.Length);

                    for (int i = 0; i < shrinkedRetainedIndices.Length; i++)
                    {
                        // Store reverse mapping to restore original dataframe indices from indices in shrinked row
                        var originalIndex = shrinkedRetainedIndices[i];
                        newRetainedIndicesReverseMapping.Add(i, originalIndex);
                    }

                    retainedIndicesReverseMapping = newRetainedIndicesReverseMapping;

                    var indices = new Int64DataFrameColumn("Indices", shrinkedRetainedIndices);
                    shrinkedRetainedColumn = shrinkedRetainedColumn.Clone(indices);
                }

                DataFrameColumn supplementaryColumn = supplementaryDataFrame.Columns[supplemetaryJoinColumnNames[colNameIndex]];

                // Find occurrenses on current step (join column)
                var newOccurrences = shrinkedRetainedColumn.GetGroupedOccurrences(supplementaryColumn, out HashSet<long> supplementaryColumnNullIndices);

                // Convert indices from in key from local (shrinked row) to indices in original dataframe
                if (retainedIndicesReverseMapping != null)
                    newOccurrences = newOccurrences.ToDictionary(kvp => retainedIndicesReverseMapping[kvp.Key], kvp => kvp.Value);

                supplementaryJoinColumnsNullIndices.UnionWith(supplementaryColumnNullIndices);

                // Shrink join result on current column by previous join columns (if any)
                // (we have to remove occurrences that doesn't exist in previous columns, because JOIN happens only if ALL left and right columns in JOIN are matched)
                if (occurrences != null)
                {
                    newOccurrences = GetShrinkedOccurences(occurrences, newOccurrences);
                }

                occurrences = newOccurrences;
            }

            return occurrences;
        }

        private static Dictionary<long, ICollection<long>> GetShrinkedOccurences(Dictionary<long, ICollection<long>> occurrences,
            Dictionary<long, ICollection<long>> newOccurrences)
        {
            var length = GetOccurencesLength(occurrences, newOccurrences);

            var shrinkedOccurences = new Dictionary<long, ICollection<long>>(length);

            foreach (var newOccurrence in newOccurrences)
            {
                var newOccurrenceKey = newOccurrence.Key;

                var list1 = (IReadOnlyList<long>)occurrences[newOccurrenceKey];
                var list2 = (IReadOnlyList<long>)newOccurrence.Value;

                var crossing = DataFrameJoinExtensions.GetSortedListsIntersection(list1, list2);

                if (crossing.Any())
                {
                    shrinkedOccurences.Add(newOccurrenceKey, crossing);
                }
            }

            return shrinkedOccurences;
        }

        private static int GetOccurencesLength(Dictionary<long, ICollection<long>> occurrences,
            Dictionary<long, ICollection<long>> newOccurrences)
        {
            var length = 0;

            foreach (var newOccurrence in newOccurrences)
            {
                var newOccurrenceKey = newOccurrence.Key;

                var list1 = (IReadOnlyList<long>)occurrences[newOccurrenceKey];
                var list2 = (IReadOnlyList<long>)newOccurrence.Value;

                var crossing = DataFrameJoinExtensions.GetSortedListsIntersection(list1, list2);

                if (crossing.Any())
                    length++;
            }

            return length;
        }

        private static HashSet<long> PerformMerging(DataFrame retainedDataFrame, string[] retainedJoinColumnNames,
            Dictionary<long, ICollection<long>> occurrences, HashSet<long> supplementaryJoinColumnsNullIndices,
            out PrimitiveDataFrameColumn<long> retainedRowIndices, out PrimitiveDataFrameColumn<long> supplementaryRowIndices,
            bool isLeftDataFrameRetained, bool isInner, bool calculateIntersection, Func<long?, long?, bool> shouldAcceptRow = null)
        {
            var retainJoinColumns = retainedJoinColumnNames.Select(name => retainedDataFrame.Columns[name]).ToList();

            var indexColumnLength = CalculateIndexColumnLength(retainedDataFrame, retainJoinColumns, occurrences,
                supplementaryJoinColumnsNullIndices, isLeftDataFrameRetained, isInner, shouldAcceptRow);

            retainedRowIndices = new Int64DataFrameColumn("RetainedIndices", indexColumnLength);
            supplementaryRowIndices = new Int64DataFrameColumn("SupplementaryIndices", indexColumnLength);

            HashSet<long> intersection = calculateIntersection ? new HashSet<long>() : null;

            var columnIndex = 0;

            for (long i = 0; i < retainedDataFrame.Columns.RowCount; i++)
            {
                long? leftDataFrameIndex = isLeftDataFrameRetained ? i : null;
                long? rightDataFrameIndex = isLeftDataFrameRetained ? null : i;

                if (!IsAnyNullValueInColumns(retainJoinColumns, i))
                {
                    // Get all row indexes from supplementary dataframe that satisfy JOIN condition
                    if (occurrences.TryGetValue(i, out ICollection<long> rowIndices))
                    {
                        foreach (long supplementaryRowIndex in rowIndices)
                        {
                            if (shouldAcceptRow?.Invoke(leftDataFrameIndex, rightDataFrameIndex) ?? true)
                            {
                                retainedRowIndices[columnIndex] = i;
                                supplementaryRowIndices[columnIndex] = supplementaryRowIndex;
                                columnIndex++;

                                // Store intersection if required
                                if (calculateIntersection)
                                {
                                    if (!intersection.Contains(supplementaryRowIndex))
                                    {
                                        intersection.Add(supplementaryRowIndex);
                                    }
                                }
                            }
                        }
                    }
                    else
                    {
                        if (isInner)
                            continue;

                        if (shouldAcceptRow?.Invoke(leftDataFrameIndex, rightDataFrameIndex) ?? true)
                        {
                            retainedRowIndices[columnIndex] = i;
                            supplementaryRowIndices[columnIndex] = null;
                            columnIndex++;
                        }
                    }
                }
                else
                {
                    foreach (long row in supplementaryJoinColumnsNullIndices)
                    {
                        if (shouldAcceptRow?.Invoke(leftDataFrameIndex, rightDataFrameIndex) ?? true)
                        {
                            retainedRowIndices[columnIndex] = i;
                            supplementaryRowIndices[columnIndex] = row;
                            columnIndex++;
                        }
                    }
                }
            }

            return intersection;
        }

        private static bool ArrangeDataFrames(DataFrame current, DataFrame other, JoinAlgorithm joinAlgorithm,
            out DataFrame retainedDataFrame, out DataFrame supplementaryDataFrame)
        {
            if (joinAlgorithm == JoinAlgorithm.Left || joinAlgorithm == JoinAlgorithm.Right)
            {
                var isLeftDataFrameRetained = joinAlgorithm == JoinAlgorithm.Left;

                retainedDataFrame = isLeftDataFrameRetained ? current : other;
                supplementaryDataFrame = isLeftDataFrameRetained ? other : current;

                return isLeftDataFrameRetained;
            }
            else if (joinAlgorithm == JoinAlgorithm.Inner)
            {
                // Use as supplementary (for Hashing) the dataframe with the smaller RowCount
                var isLeftDataFrameRetained = current.Rows.Count > other.Rows.Count;

                retainedDataFrame = isLeftDataFrameRetained ? current : other;
                supplementaryDataFrame = isLeftDataFrameRetained ? other : current;

                return isLeftDataFrameRetained;
            }
            else if (joinAlgorithm == JoinAlgorithm.FullOuter)
            {
                // In full outer join we would like to retain data from both sides,
                // So we do it into 2 steps: one first we do LEFT JOIN and then add lost data from the RIGHT side

                retainedDataFrame = current;
                supplementaryDataFrame = other;

                return true;
            }
            else
            {
                throw new NotImplementedException(nameof(joinAlgorithm));
            }
        }

        private static int CalculateIndexColumnLength(DataFrame retainedDataFrame, IReadOnlyCollection<DataFrameColumn> retainJoinColumns,
            Dictionary<long, ICollection<long>> occurrences, HashSet<long> supplementaryJoinColumnsNullIndices,
            bool isLeftDataFrameRetained, bool isInner, Func<long?, long?, bool> shouldAcceptRow = null)
        {
            var indexColumnLength = 0;

            for (long i = 0; i < retainedDataFrame.Columns.RowCount; i++)
            {
                long? leftDataFrameIndex = isLeftDataFrameRetained ? i : null;
                long? rightDataFrameIndex = isLeftDataFrameRetained ? null : i;

                if (!IsAnyNullValueInColumns(retainJoinColumns, i))
                {
                    // Get all row indexes from supplementary dataframe that satisfy JOIN condition
                    if (occurrences.TryGetValue(i, out ICollection<long> rowIndices))
                    {
                        if (shouldAcceptRow == null)
                        {
                            indexColumnLength += rowIndices.Count;
                        }
                        else
                        {
                            foreach (long supplementaryRowIndex in rowIndices)
                            {
                                if (shouldAcceptRow(leftDataFrameIndex, rightDataFrameIndex))
                                {
                                    indexColumnLength++;
                                }
                            }
                        }
                    }
                    else
                    {
                        if (isInner)
                            continue;

                        if (shouldAcceptRow?.Invoke(leftDataFrameIndex, rightDataFrameIndex) ?? true)
                        {
                            indexColumnLength++;
                        }
                    }
                }
                else
                {
                    if (shouldAcceptRow == null)
                    {
                        indexColumnLength += supplementaryJoinColumnsNullIndices.Count;
                    }
                    else
                    {
                        foreach (long row in supplementaryJoinColumnsNullIndices)
                        {
                            if (shouldAcceptRow(leftDataFrameIndex, rightDataFrameIndex))
                            {
                                indexColumnLength++;
                            }
                        }
                    }
                }
            }

            return indexColumnLength;
        }

        public DataFrame Merge(DataFrame other, string[] leftJoinColumns, string[] rightJoinColumns, string leftSuffix = "_left", string rightSuffix = "_right",
            JoinAlgorithm joinAlgorithm = JoinAlgorithm.Left, IReadOnlyCollection<string> leaveColumns = null,
            Func<long?, long?, bool> shouldAcceptRow = null)
        {
            if (other == null)
                throw new ArgumentNullException(nameof(other));

            // In Outer join the joined dataframe retains each row — even if no other matching row exists in supplementary dataframe.
            // Outer joins subdivide further into left outer joins (left dataframe is retained), right outer joins (rightdataframe is retained), in full outer both are retained

            DataFrame retainedDataFrame;
            DataFrame supplementaryDataFrame;
            PrimitiveDataFrameColumn<long> retainedRowIndices;
            PrimitiveDataFrameColumn<long> supplementaryRowIndices;

            bool isLeftDataFrameRetained = ArrangeDataFrames(this, other, joinAlgorithm,
                out retainedDataFrame, out supplementaryDataFrame);

            var supplementaryJoinColumns = isLeftDataFrameRetained ? rightJoinColumns : leftJoinColumns;
            var retainedJoinColumns = isLeftDataFrameRetained ? leftJoinColumns : rightJoinColumns;

            var intersection = Merge(retainedDataFrame, supplementaryDataFrame, retainedJoinColumns, supplementaryJoinColumns,
                    out retainedRowIndices, out supplementaryRowIndices, isLeftDataFrameRetained,
                    isInner: joinAlgorithm == JoinAlgorithm.Inner,
                    calculateIntersection: joinAlgorithm == JoinAlgorithm.FullOuter,
                    shouldAcceptRow: shouldAcceptRow);

            if (joinAlgorithm == JoinAlgorithm.FullOuter)
            {
                // Step 1
                // LEFT JOIN have done in merge and got intersection

                // Step 2
                // Do RIGHT JOIN to retain all data from supplementary DataFrame too (take into account data intersection from the first step to avoid duplicates)

                var columns = supplementaryJoinColumns.Select(name => supplementaryDataFrame.Columns[name]).ToList();

                for (long i = 0; i < supplementaryDataFrame.Columns.RowCount; i++)
                {
                    if (!IsAnyNullValueInColumns(columns, i))
                    {
                        if (!intersection.Contains(i))
                        {
                            retainedRowIndices.Append(null);
                            supplementaryRowIndices.Append(i);
                        }
                    }
                }
            }

            DataFrame ret = new DataFrame();

            PrimitiveDataFrameColumn<long> mapIndicesLeft = isLeftDataFrameRetained ? retainedRowIndices : supplementaryRowIndices;
            PrimitiveDataFrameColumn<long> mapIndicesRight = isLeftDataFrameRetained ? supplementaryRowIndices : retainedRowIndices;

            if (leaveColumns?.Any() ?? false)
            {
                var retColumnIndex = 0;

                foreach (var column in this.Columns.Where(x => leaveColumns.Contains(x.Name)))
                {
                    var clone = column.Clone(mapIndicesLeft);
                    clone.SetName(clone.Name + leftSuffix);

                    ret.Columns.Insert(retColumnIndex, clone);
                    retColumnIndex++;
                }

                foreach (var column in other.Columns.Where(x => leaveColumns.Contains(x.Name)))
                {
                    var clone = column.Clone(mapIndicesRight);
                    clone.SetName(clone.Name + rightSuffix);

                    ret.Columns.Insert(retColumnIndex, clone);
                    retColumnIndex++;
                }
            }
            else
            {
                // Insert columns from left dataframe (this)
                for (int i = 0; i < this.Columns.Count; i++)
                {
                    ret.Columns.Insert(i, this.Columns[i].Clone(mapIndicesLeft));
                }

                // Insert columns from right dataframe (other)
                for (int i = 0; i < other.Columns.Count; i++)
                {
                    DataFrameColumn column = other.Columns[i].Clone(mapIndicesRight);

                    SetSuffixForDuplicatedColumnNames(ret, column, leftSuffix, rightSuffix);
                    ret.Columns.Insert(ret.Columns.Count, column);
                }
            }

            return ret;
        }
    }
}
