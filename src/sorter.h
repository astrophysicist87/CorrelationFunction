#ifndef SORTER_H
#define SORTER_H

#include <vector>
#include <algorithm>

#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>

#define COLUMN_INDEX_TO_SORT_BY	1	// duh

using namespace std;

template <typename T>
vector<size_t> ordered(vector<T> const& values)
    {
        using namespace boost::phoenix;
        using namespace boost::phoenix::arg_names;

        vector<size_t> indices(values.size());
        int i = 0;
        transform(values.begin(), values.end(), indices.begin(), ref(i)++);
        sort(indices.begin(), indices.end(), ref(values)[arg1] < ref(values)[arg2]);
        return indices;
    }

template <typename T>
bool cmp(const vector<T>& a, const vector<T>& b)
{
    return a[COLUMN_INDEX_TO_SORT_BY] < b[COLUMN_INDEX_TO_SORT_BY];
}

template <typename T>
void sort_by_column(vector< vector<T> > * values)
{
	std::stable_sort((*values).begin(), (*values).end(), cmp<T>);

	return;
}

void inline debugger(int cln, const char* cfn)
{
	cerr << "You made it to " << cfn << ":" << cln << "!" << endl;
	return;
}

#endif
