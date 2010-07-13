/* AUTHOR: Vladislav S. Yakovlev */

#include "reader.hh"
#include <cstdlib>
#include <string>

//------------------------------------------------------------------

void ReadTable(istream& is, vector< vector<Real> >& table)
{
    int i, n;
//    int line_counter;
    const int max_string_length = 1024;
    char str1[max_string_length];
    istringstream str_stream;
    Real x;
    vector<Real> row;

    n = table.size();
    for (i=0; i<n; i++) table[i].clear();
    table.clear();

    // read from the input stream line by line
//    line_counter = 0;
    while (is.getline(str1, max_string_length))
    {
//	line_counter++;
	row.clear();
	// if the first non-blank character is '#', it's a comment string
	// we also should check whether it's an empty string
	for (i=0; i<max_string_length && str1[i]!='\0' && str1[i]!='\n' &&
		 (str1[i]==' ' || str1[i]=='\t'); i++) {}
	if (i==max_string_length) continue;
	if (str1[i]=='\0' || str1[i]=='\n' || str1[i]=='#') continue;
	// try to interpret the string
	str_stream.str(string(str1));
	str_stream.clear();
	while(str_stream >> x) row.push_back(x);
	if (row.size() == 0)
	{
	    string message(str1);
	    message = "\
can't interpret the following line as whitespace-separated real numbers:\n" +
		message;
	    throw message;
	}
	table.push_back(row);
    }
}

//------------------------------------------------------------------
