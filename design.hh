#ifndef __DESIGN_HH
#define __DESIGN_HH

#include "common.hh"
#include <string>
#include <vector>

/* An object of the class 'Design' knows how to read itself from a file, how
   to write itself to a file, and it can tell some simple information about
   the design. An object of the class 'Coating' can do more.
*/
class Design
{
public:
    Design();
    ~Design();
    string name;
    void ReadFromFile(const string& file_name);

    void SaveToFile(string file_name="");
    // the name is not given, it is taken from the member string 'name'

    int NumberOfLayers() const {return no_of_layers;}
    bool isLike(const Design&);
    /* returns true if the given design has the same number of layers made
       from exactly the same materials; of some of the layers are "fixed",
       they must be fixed in both designs, and their thicknesses must be
       the same */
    void Clear(); // removes all the layers
    void AddLayer(const string& material_name, Real thickness,
		  bool fixed_thickness=false); // on the top
    const string& MaterialName(int layer_index) const;
    Real LayerThickness(int layer_index) const;
    bool FixedThickness(int layer_index) const;

private:
    vector<string> material_names;
    vector<Real> thicknesses;
    vector<bool> fixed_thicknesses; // marks the layers, which may not
				    // be optimised
    int no_of_layers;
};

//-------------------------------------------------------------------------

void RenameDesignsToCWD(vector<Design>& designs);
/* extracts the path name from each design name (if present) and makes sure
   that no two designs get the same name; if they do, the function changes
   the conflicting names */


#endif
