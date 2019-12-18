#include "testtecplot.h"

#ifdef DEBUG
template <> int clVector_namespace::clVector<int>::_num_vector=0;
template <> int clVector_namespace::clVector<double>::_num_vector=0;
template <> int clVector_namespace::clMatrix<int>::_num_matrix=0;
template <> int clVector_namespace::clMatrix<double>::_num_matrix=0;
template <> int clVector_namespace::clMatrix3D<double>::_num_matrix=0;
int clString_namespace::clString::_num_string=0;
int clString_namespace::clString_vector::_num_string_vector=0;
#endif

const int test_bin=1;
const int test_write=1;
const int test_eliminate_varsharing=1;
const int test_find_minmax_vertex=1;
const int test_translate=1;
const int test_add_var_and_copy=1;
const int test_add_data=1;
const int test_add_remove_var=1;
const int test_combine_zones=0;
const int test_interpolation_1=1;
const int test_interpolation_2=1;
const int test_interpolation_3=1;
const int test_interpolation_4=1;
const int test_interpolation_5=1;
const int test_interpolation_6=1;
const int test_interpolation_cc=1;
const int test_copy_constructor_1=1;
const int test_copy_constructor_2=1;
const int test_copy_constructor_3=1;

int main(int argc, char* argv[])
{
	const char inputdir[]="Input";
	const char outputdir[]="Output";
	char filename[NAME_LENGTH];
	clTecplot *tecplot=NULL,*ttecplot=NULL;

	/* Test read structured */
	if (test_bin) {
		tecplot=new clTecplot;
		sprintf(filename,"%s/test_structured_1zone.plt",inputdir);
		tecplot->read(filename);
		sprintf(filename,"%s/output_test_bin_structured_1zone.dat",outputdir);
		tecplot->write_ascii(filename);
		sprintf(filename,"%s/output_test_bin_structured_1zone.plt",outputdir);
		tecplot->write(filename);
		delete tecplot;
		tecplot=new clTecplot;
		sprintf(filename,"%s/test_structured_2zones.plt",inputdir);
		tecplot->read(filename);
		sprintf(filename,"%s/output_test_bin_structured_2zones.dat",outputdir);
		tecplot->write_ascii(filename);
		sprintf(filename,"%s/output_test_bin_structured_2zones.plt",outputdir);
		tecplot->write(filename);
		delete tecplot;
		tecplot=new clTecplot;
		sprintf(filename,"%s/test_structured_block.plt",inputdir);
		tecplot->read(filename);
		sprintf(filename,"%s/output_test_bin_structured_block.dat",outputdir);
		tecplot->write_ascii(filename);
		sprintf(filename,"%s/output_test_bin_structured_block.plt",outputdir);
		tecplot->write(filename);
		delete tecplot;
	}
	/* Test read unstructured */
	if (test_bin) {
		tecplot=new clTecplot;
		sprintf(filename,"%s/test_unstructured_1zone.plt",inputdir);
		tecplot->read(filename);
		sprintf(filename,"%s/output_test_bin_unstructured_1zone.dat",outputdir);
		tecplot->write_ascii(filename);
		sprintf(filename,"%s/output_test_bin_unstructured_1zone.plt",outputdir);
		tecplot->write(filename);
		delete tecplot;
		tecplot=new clTecplot;
		sprintf(filename,"%s/test_unstructured_2zones.plt",inputdir);
		tecplot->read(filename);
		sprintf(filename,"%s/output_test_bin_unstructured_2zones.dat",outputdir);
		tecplot->write_ascii(filename);
		sprintf(filename,"%s/output_test_bin_unstructured_2zones.plt",outputdir);
		tecplot->write(filename);
		delete tecplot;
		tecplot=new clTecplot;
		sprintf(filename,"%s/test_unstructured_2zones_cc.plt",inputdir);
		tecplot->read(filename);
		sprintf(filename,"%s/output_test_bin_unstructured_2zones_cc.dat",outputdir);
		tecplot->write_ascii(filename);
		sprintf(filename,"%s/output_test_bin_unstructured_2zones_cc.plt",outputdir);
		tecplot->write(filename);
		delete tecplot;
		tecplot=new clTecplot;
		sprintf(filename,"%s/test_unstructured_cc.plt",inputdir);
		tecplot->read(filename);
		sprintf(filename,"%s/output_test_bin_unstructured_cc.dat",outputdir);
		tecplot->write_ascii(filename);
		sprintf(filename,"%s/output_test_bin_unstructured_cc.plt",outputdir);
		tecplot->write(filename);
		delete tecplot;
		tecplot=new clTecplot;
		sprintf(filename,"%s/test_unstructured_shared_quad.plt",inputdir);
		tecplot->read(filename);
		sprintf(filename,"%s/output_test_bin_unstructured_shared_quad.dat",outputdir);
		tecplot->write_ascii(filename);
		sprintf(filename,"%s/output_test_bin_unstructured_shared_quad.plt",outputdir);
		tecplot->write(filename);
		delete tecplot;
	}

	/* Test write_ascii structured */
	if (test_write) {
		tecplot=new clTecplot;
		sprintf(filename,"%s/test_structured_1zone.dat",inputdir);
		tecplot->read(filename);
		sprintf(filename,"%s/output_test_write_structured_1zone.dat",outputdir);
		tecplot->write_ascii(filename);
		delete tecplot;
		tecplot=new clTecplot;
		sprintf(filename,"%s/test_structured_2zones.dat",inputdir);
		tecplot->read(filename);
		sprintf(filename,"%s/output_test_write_structured_2zones.dat",outputdir);
		tecplot->write_ascii(filename);
		delete tecplot;
		tecplot=new clTecplot;
		sprintf(filename,"%s/test_structured_2zones.dat",inputdir);
		tecplot->read(filename);
		sprintf(filename,"%s/output_test_write_structured_2zones_block.dat",outputdir);
		tecplot->write_ascii(filename,"BLOCK");
		delete tecplot;
		tecplot=new clTecplot;
		sprintf(filename,"%s/test_structured_noijk.dat",inputdir);
		tecplot->read(filename);
		sprintf(filename,"%s/output_test_write_structured_noijk.dat",outputdir);
		tecplot->write_ascii(filename);
		delete tecplot;
		tecplot=new clTecplot;
	}

	/* Test write_ascii unstructured */
	if (test_write) {
		tecplot=new clTecplot;
		sprintf(filename,"%s/test_unstructured_1zone.dat",inputdir);
		tecplot->read(filename);
		sprintf(filename,"%s/output_test_write_unstructured_1zone.dat",outputdir);
		tecplot->write_ascii(filename);
		delete tecplot;
		tecplot=new clTecplot;
		sprintf(filename,"%s/test_unstructured_2zones.dat",inputdir);
		tecplot->read(filename);
		sprintf(filename,"%s/output_test_write_unstructured_2zones.dat",outputdir);
		tecplot->write_ascii(filename);
		delete tecplot;
		tecplot=new clTecplot;
		sprintf(filename,"%s/test_unstructured_2zones.dat",inputdir);
		tecplot->read(filename);
		sprintf(filename,"%s/output_test_write_unstructured_2zones_block.dat",outputdir);
		tecplot->write_ascii(filename,"BLOCK");
		delete tecplot;
		tecplot=new clTecplot;
		sprintf(filename,"%s/test_unstructured_cc.dat",inputdir);
		tecplot->read(filename);
		sprintf(filename,"%s/output_test_write_unstructured_cc.dat",outputdir);
		tecplot->write_ascii(filename,"BLOCK");
		delete tecplot;
		tecplot=new clTecplot;
		sprintf(filename,"%s/test_unstructured_shared_quad.dat",inputdir);
		tecplot->read(filename);
		sprintf(filename,"%s/output_test_write_unstructured_shared_quad.dat",outputdir);
		tecplot->write_ascii(filename,"BLOCK");
		delete tecplot;
	}

	/* Test eliminate varsharing and convert_fequadrilateral_to_fetriangle */
	if (test_eliminate_varsharing) {
		tecplot=new clTecplot;
		sprintf(filename,"%s/test_unstructured_shared_quad.dat",inputdir);
		tecplot->read(filename);
		tecplot->convert_fequadrilateral_to_fetriangle();
		tecplot->eliminate_varsharing();
		sprintf(filename,"%s/output_test_eliminate_varsharing_unstructured_unshared_tri.dat",outputdir);
		tecplot->write_ascii(filename);
		delete tecplot;
	}

	/* Test translate data */
	if (test_translate) {
		tecplot=new clTecplot;
		sprintf(filename,"%s/test_structured_1zone.dat",inputdir);
		tecplot->read(filename);
		tecplot->translate_data(2.0,0);
		tecplot->translate_data(-2.0,1);
		sprintf(filename,"%s/output_test_translate_structured_1zone.dat",outputdir);
		tecplot->write_ascii(filename);
		delete tecplot;
		tecplot=new clTecplot;
		sprintf(filename,"%s/test_unstructured_1zone.dat",inputdir);
		tecplot->read(filename);
		tecplot->translate_data(2.0,0);
		tecplot->translate_data(-2.0,1);
		sprintf(filename,"%s/output_test_translate_unstructured_1zone.dat",outputdir);
		tecplot->write_ascii(filename);
		delete tecplot;
	}

	/* Test add variable and copy */
	if (test_add_var_and_copy) {
		tecplot=new clTecplot;
		sprintf(filename,"%s/test_structured_2zones.dat",inputdir);
		tecplot->read(filename);
		ttecplot=new clTecplot;
		ttecplot->copy(tecplot);
		sprintf(filename,"%s/output_test_copy_structured_2zones.dat",outputdir);
		ttecplot->write_ascii(filename);
		ttecplot->add_variable();
		ttecplot->copy_variable(tecplot,4,5);
		sprintf(filename,"%s/output_test_add_and_copy_var_structured_2zones.dat",outputdir);
		ttecplot->write_ascii(filename);
		delete tecplot;
		delete ttecplot;
		tecplot=new clTecplot;
		sprintf(filename,"%s/test_unstructured_2zones.dat",inputdir);
		tecplot->read(filename);
		ttecplot=new clTecplot;
		ttecplot->copy(tecplot);
		sprintf(filename,"%s/output_test_copy_unstructured_2zones.dat",outputdir);
		ttecplot->write_ascii(filename);
		ttecplot->add_variable();
		ttecplot->copy_variable(tecplot,3,4);
		sprintf(filename,"%s/output_test_add_and_copy_var_unstructured_2zones.dat",outputdir);
		ttecplot->write_ascii(filename);
		delete tecplot;
		delete ttecplot;
	}

	/* Test add data */
	if (test_add_data) {
		tecplot=new clTecplot;
		sprintf(filename,"%s/test_structured_2zones.dat",inputdir);
		tecplot->read(filename);
		ttecplot=new clTecplot;
		ttecplot->copy(tecplot);
		tecplot->add_data(ttecplot);
		sprintf(filename,"%s/output_test_add_data_structured_2zones.dat",outputdir);
		tecplot->write_ascii(filename);
		delete tecplot;
		delete ttecplot;
		tecplot=new clTecplot;
		sprintf(filename,"%s/test_unstructured_2zones.dat",inputdir);
		tecplot->read(filename);
		ttecplot=new clTecplot;
		ttecplot->copy(tecplot);
		tecplot->add_data(ttecplot);
		sprintf(filename,"%s/output_test_add_data_unstructured_2zones.dat",outputdir);
		tecplot->write_ascii(filename);
		delete tecplot;
		delete ttecplot;
	}

	/* Test add/remove variable */
	if (test_add_remove_var) {
		tecplot=new clTecplot;
		sprintf(filename,"%s/test_structured_2zones.dat",inputdir);
		tecplot->read(filename);
		tecplot->remove_variable(3);
		tecplot->remove_variable(1);
		tecplot->add_variables(3);
		sprintf(filename,"%s/output_test_add_remove_var_structured_2zones.dat",outputdir);
		tecplot->write_ascii(filename);
		delete tecplot;
		tecplot=new clTecplot;
		sprintf(filename,"%s/test_unstructured_2zones.dat",inputdir);
		tecplot->read(filename);
		tecplot->remove_variable(3);
		tecplot->remove_variable(1);
		tecplot->add_variables(3);
		sprintf(filename,"%s/output_test_add_remove_var_unstructured_2zones.dat",outputdir);
		tecplot->write_ascii(filename);
		delete tecplot;
	}

	/* Test combine zones: Make this into an internal clTecplot routine */
	if (test_combine_zones) {
		tecplot=new clTecplot;
		sprintf(filename,"%s/test_unstructured_2zones_cc.dat",inputdir);
		tecplot->read(filename);
		int n_data,num_data,num_elem=0,num_node=0,n_var,n_zone;
		int num_zone=tecplot->get_num_zone();
		const int num_var=tecplot->get_num_var();
		clTecplot_zone *currentzone=NULL;
		for (n_zone=0; n_zone<num_zone; n_zone++) {
			currentzone=tecplot->get_zoneptr(n_zone);
			num_node+=currentzone->get_num_node();
			num_elem+=currentzone->get_num_elem();
		}
		clVector_namespace::clVector<double> *data=new clVector_namespace::clVector<double> [num_var];
		double **dataptr=new double* [num_var];
		for (n_var=0; n_var<num_var; n_var++) {
			currentzone=tecplot->get_zoneptr(0);
			if (!currentzone->get_varlocation(n_var)) num_data=num_node;
			else num_data=num_elem;
			data[n_var].allocate(num_data);
			dataptr[n_var]=data[n_var].get_startptr();
		}
		const int num_conn=(tecplot->get_zoneptr(0))->get_num_conn();
		clVector_namespace::clMatrix<int> conn(num_elem,num_conn);
		int num_node_offset=0;
		int *connptr=conn.get_startptr();
		const int *connorgptr=NULL;
		const double *dataorgptr=NULL;
		for (n_zone=0; n_zone<num_zone; n_zone++) {
			for (n_var=0; n_var<num_var; n_var++) {
				if (!currentzone->get_varlocation(n_var)) num_data=currentzone->get_num_node();
				else num_data=currentzone->get_num_elem();
				currentzone=tecplot->get_zoneptr(n_zone);
				dataorgptr=currentzone->get_datastartptr(n_var);
				for (n_data=0; n_data<num_data; n_data++) {
					*dataptr[n_var]=(*dataorgptr++);
					dataptr[n_var]++;
				}
			}
			connorgptr=currentzone->get_connstartptr();
			num_data=num_conn*currentzone->get_num_elem();
			for (n_data=0; n_data<num_data; n_data++) (*connptr++)=(*connorgptr++)+num_node_offset;
			num_node_offset+=currentzone->get_num_node();
		}
		tecplot->remove_zone(1,num_zone-1);
		currentzone=tecplot->get_zoneptr(0);
		currentzone->set_data(data);
		currentzone->set_conn(conn);
		currentzone->set_num_node(num_node);
		currentzone->set_num_elem(num_elem);
		currentzone->set_title("test_combined");
		sprintf(filename,"%s/output_test_combine_zones_unstructured_cc_2zones.dat",outputdir);
		tecplot->write_ascii(filename);
		delete [] data;
		delete [] dataptr;
	}

	/* Test interpolation */
	if (test_interpolation_1) {
		tecplot=new clTecplot;
		ttecplot=new clTecplot;
		sprintf(filename,"%s/test_unstructured_1zone.dat",inputdir);
		tecplot->read(filename);
		sprintf(filename,"%s/test_unstructured_1zone_fine.dat",inputdir);
		ttecplot->read(filename);
		ttecplot->remove_variable(4);
		ttecplot->remove_variable(3);
		ttecplot->add_variables(2);
		const int num_zone=tecplot->get_num_zone();
		if (num_zone!=ttecplot->get_num_zone()) clUtils_namespace::fatal_err("num_zone!=tecplot->get_num_zone()");
		const clTecplot_zone *tecplotzone=NULL;
		clTecplot_zone *ttecplotzone=NULL;
		for (int n_zone=0; n_zone<num_zone; n_zone++) {
			tecplotzone=tecplot->get_zoneptr(n_zone);
			ttecplotzone=ttecplot->get_zoneptr(n_zone);
			tecplotzone->interpolate(ttecplotzone);
		}
		sprintf(filename,"%s/output_test_interpolation_1_unstructured_1zone.dat",outputdir);
		ttecplot->write_ascii(filename);
		delete tecplot;
		delete ttecplot;
	}
	if (test_interpolation_2) {
		tecplot=new clTecplot;
		ttecplot=new clTecplot;
		sprintf(filename,"%s/test_unstructured_noconn.dat",inputdir);
		tecplot->read(filename);
		sprintf(filename,"%s/test_unstructured_noconn_fine.dat",inputdir);
		ttecplot->read(filename);
		ttecplot->remove_variable(4);
		ttecplot->remove_variable(3);
		ttecplot->add_variables(2);
		const int num_zone=tecplot->get_num_zone();
		if (num_zone!=ttecplot->get_num_zone()) clUtils_namespace::fatal_err("num_zone!=tecplot->get_num_zone()");
		const clTecplot_zone *tecplotzone=NULL;
		clTecplot_zone *ttecplotzone=NULL;
		for (int n_zone=0; n_zone<num_zone; n_zone++) {
			tecplotzone=tecplot->get_zoneptr(n_zone);
			ttecplotzone=ttecplot->get_zoneptr(n_zone);
			tecplotzone->interpolate(ttecplotzone);
		}
		sprintf(filename,"%s/output_test_interpolation_2_unstructured_noconn.dat",outputdir);
		ttecplot->write_ascii(filename);
		delete tecplot;
		delete ttecplot;
	}
	if (test_interpolation_3) {
		tecplot=new clTecplot;
		ttecplot=new clTecplot;
		sprintf(filename,"%s/test_unstructured_1zone.dat",inputdir);
		tecplot->read(filename);
		ttecplot->read(filename);
		ttecplot->remove_variable(4);
		ttecplot->remove_variable(3);
		ttecplot->add_variables(1);
		ttecplot->set_varname(3,tecplot->get_varname(4));
		const clVector_namespace::clVector<int> usevar(4,-1);
		usevar(3)=4;
		const int num_zone=tecplot->get_num_zone();
		if (num_zone!=ttecplot->get_num_zone()) clUtils_namespace::fatal_err("num_zone!=tecplot->get_num_zone()");
		const clTecplot_zone *tecplotzone=NULL;
		clTecplot_zone *ttecplotzone=NULL;
		for (int n_zone=0; n_zone<num_zone; n_zone++) {
			tecplotzone=tecplot->get_zoneptr(n_zone);
			ttecplotzone=ttecplot->get_zoneptr(n_zone);
			tecplotzone->interpolate(ttecplotzone,usevar);
		}
		sprintf(filename,"%s/output_test_interpolation_3_unstructured_1zone.dat",outputdir);
		ttecplot->write_ascii(filename);
		delete tecplot;
		delete ttecplot;
	}
	if (test_interpolation_4) {
		tecplot=new clTecplot;
		ttecplot=new clTecplot;
		sprintf(filename,"%s/test_unstructured_1zone.dat",inputdir);
		tecplot->read(filename);
		tecplot->remove_variable(3);
		ttecplot->read(filename);
		ttecplot->remove_variable(4);
		ttecplot->add_variables(1);
		ttecplot->set_varname(4,tecplot->get_varname(3));
		const clVector_namespace::clVector<int> usevar(5,-1);
		usevar(4)=3;
		const int num_zone=tecplot->get_num_zone();
		if (num_zone!=ttecplot->get_num_zone()) clUtils_namespace::fatal_err("num_zone!=tecplot->get_num_zone()");
		const clTecplot_zone *tecplotzone=NULL;
		clTecplot_zone *ttecplotzone=NULL;
		for (int n_zone=0; n_zone<num_zone; n_zone++) {
			tecplotzone=tecplot->get_zoneptr(n_zone);
			ttecplotzone=ttecplot->get_zoneptr(n_zone);
			tecplotzone->interpolate(ttecplotzone,usevar);
		}
		sprintf(filename,"%s/output_test_interpolation_4_unstructured_1zone.dat",outputdir);
		ttecplot->write_ascii(filename);
		delete tecplot;
		delete ttecplot;
	}
	if (test_interpolation_5) {
		tecplot=new clTecplot;
		ttecplot=new clTecplot;
		sprintf(filename,"%s/test_unstructured_cc_fine.dat",inputdir);
		tecplot->read(filename);
		ttecplot->read(filename);
		ttecplot->remove_variable(4);
		ttecplot->remove_variable(3);
		ttecplot->add_variables(1,1);
		ttecplot->set_varname(3,tecplot->get_varname(4));
		const clVector_namespace::clVector<int> usevar(4,-1);
		usevar(3)=4;
		const int num_zone=tecplot->get_num_zone();
		if (num_zone!=ttecplot->get_num_zone()) clUtils_namespace::fatal_err("num_zone!=tecplot->get_num_zone()");
		const clTecplot_zone *tecplotzone=NULL;
		clTecplot_zone *ttecplotzone=NULL;
		for (int n_zone=0; n_zone<num_zone; n_zone++) {
			tecplotzone=tecplot->get_zoneptr(n_zone);
			ttecplotzone=ttecplot->get_zoneptr(n_zone);
			tecplotzone->interpolate(ttecplotzone,usevar);
		}
		sprintf(filename,"%s/output_test_interpolation_5_unstructured_cc_fine.dat",outputdir);
		ttecplot->write_ascii(filename);
		delete tecplot;
		delete ttecplot;
	}
	if (test_interpolation_6) {
		tecplot=new clTecplot;
		ttecplot=new clTecplot;
		sprintf(filename,"%s/test_unstructured_cc_fine.dat",inputdir);
		tecplot->read(filename);
		tecplot->remove_variable(3);
		ttecplot->read(filename);
		ttecplot->remove_variable(4);
		ttecplot->add_variables(1,1);
		ttecplot->set_varname(4,tecplot->get_varname(3));
		const clVector_namespace::clVector<int> usevar(5,-1);
		usevar(4)=3;
		const int num_zone=tecplot->get_num_zone();
		if (num_zone!=ttecplot->get_num_zone()) clUtils_namespace::fatal_err("num_zone!=tecplot->get_num_zone()");
		const clTecplot_zone *tecplotzone=NULL;
		clTecplot_zone *ttecplotzone=NULL;
		for (int n_zone=0; n_zone<num_zone; n_zone++) {
			tecplotzone=tecplot->get_zoneptr(n_zone);
			ttecplotzone=ttecplot->get_zoneptr(n_zone);
			tecplotzone->interpolate(ttecplotzone,usevar);
		}
		sprintf(filename,"%s/output_test_interpolation_6_unstructured_cc_fine.dat",outputdir);
		ttecplot->write_ascii(filename);
		delete tecplot;
		delete ttecplot;
	}

	if (test_interpolation_cc) {
		tecplot=new clTecplot;
		ttecplot=new clTecplot;
		sprintf(filename,"%s/test_unstructured_cc.dat",inputdir);
		tecplot->read(filename);
		tecplot->remove_variable(3);
		ttecplot->read(filename);
		ttecplot->remove_variable(4);
		ttecplot->add_variables(1,1);
		ttecplot->set_varname(4,tecplot->get_varname(3));
		const clVector_namespace::clVector<int> usevar(5,-1);
		usevar(4)=3;
		const int num_zone=tecplot->get_num_zone();
		if (num_zone!=ttecplot->get_num_zone()) clUtils_namespace::fatal_err("num_zone!=tecplot->get_num_zone()");
		const clTecplot_zone *tecplotzone=NULL;
		clTecplot_zone *ttecplotzone=NULL;
		for (int n_zone=0; n_zone<num_zone; n_zone++) {
			tecplotzone=tecplot->get_zoneptr(n_zone);
			ttecplotzone=ttecplot->get_zoneptr(n_zone);
			tecplotzone->interpolate_cc(ttecplotzone,usevar);
		}
		sprintf(filename,"%s/output_test_interpolation_cc_unstructured_cc.dat",outputdir);
		ttecplot->write_ascii(filename);
		delete tecplot;
		delete ttecplot;
	}

	/* Test copy constructor structured */
	if (test_copy_constructor_1) {
		tecplot=new clTecplot;
		sprintf(filename,"%s/test_structured_2zones.dat",inputdir);
		tecplot->read(filename);
		ttecplot=new clTecplot(*tecplot);
		delete tecplot;
		sprintf(filename,"%s/output_test_copy_constructor_1_structured_2zones.dat",outputdir);
		ttecplot->write_ascii(filename);
		delete ttecplot;
		tecplot=new clTecplot;
		sprintf(filename,"%s/test_unstructured_2zones.dat",inputdir);
		tecplot->read(filename);
		ttecplot=new clTecplot(*tecplot);
		delete tecplot;
		sprintf(filename,"%s/output_test_copy_constructor_1_unstructured_2zones.dat",outputdir);
		ttecplot->write_ascii(filename);
		delete ttecplot;
		tecplot=new clTecplot;
		sprintf(filename,"%s/test_unstructured_2zones_cc.dat",inputdir);
		tecplot->read(filename);
		ttecplot=new clTecplot(*tecplot);
		delete tecplot;
		sprintf(filename,"%s/output_test_copy_constructor_1_unstructured_2zones_cc.dat",outputdir);
		ttecplot->write_ascii(filename);
		delete ttecplot;
		tecplot=new clTecplot;
		sprintf(filename,"%s/test_unstructured_noconn.dat",inputdir);
		tecplot->read(filename);
		ttecplot=new clTecplot(*tecplot);
		delete tecplot;
		sprintf(filename,"%s/output_test_copy_constructor_1_unstructured_noconn.dat",outputdir);
		ttecplot->write_ascii(filename);
		delete ttecplot;
	}
	if (test_copy_constructor_2) {
		tecplot=new clTecplot;
		sprintf(filename,"%s/test_structured_2zones.dat",inputdir);
		tecplot->read(filename);
		ttecplot=new clTecplot;
		ttecplot->copy(tecplot);
		delete tecplot;
		sprintf(filename,"%s/output_test_copy_constructor_2_structured_2zones.dat",outputdir);
		ttecplot->write_ascii(filename);
		delete ttecplot;
		tecplot=new clTecplot;
		sprintf(filename,"%s/test_unstructured_2zones.dat",inputdir);
		tecplot->read(filename);
		ttecplot=new clTecplot;
		ttecplot->copy(tecplot);
		delete tecplot;
		sprintf(filename,"%s/output_test_copy_constructor_2_unstructured_2zones.dat",outputdir);
		ttecplot->write_ascii(filename);
		delete ttecplot;
		tecplot=new clTecplot;
		sprintf(filename,"%s/test_unstructured_2zones_cc.dat",inputdir);
		tecplot->read(filename);
		ttecplot=new clTecplot;
		ttecplot->copy(tecplot);
		delete tecplot;
		sprintf(filename,"%s/output_test_copy_constructor_2_unstructured_2zones_cc.dat",outputdir);
		ttecplot->write_ascii(filename);
		delete ttecplot;
		tecplot=new clTecplot;
		sprintf(filename,"%s/test_unstructured_noconn.dat",inputdir);
		tecplot->read(filename);
		ttecplot=new clTecplot(*tecplot);
		ttecplot->copy(tecplot);
		delete tecplot;
		sprintf(filename,"%s/output_test_copy_constructor_2_unstructured_noconn.dat",outputdir);
		ttecplot->write_ascii(filename);
		delete ttecplot;
	}
	if (test_copy_constructor_3) {
		tecplot=new clTecplot;
		sprintf(filename,"%s/test_structured_2zones.dat",inputdir);
		tecplot->read(filename);
		ttecplot=new clTecplot;
		(*ttecplot)=(*tecplot);
		delete tecplot;
		sprintf(filename,"%s/output_test_copy_constructor_3_structured_2zones.dat",outputdir);
		ttecplot->write_ascii(filename);
		delete ttecplot;
		tecplot=new clTecplot;
		sprintf(filename,"%s/test_unstructured_2zones.dat",inputdir);
		tecplot->read(filename);
		ttecplot=new clTecplot;
		(*ttecplot)=(*tecplot);
		delete tecplot;
		sprintf(filename,"%s/output_test_copy_constructor_3_unstructured_2zones.dat",outputdir);
		ttecplot->write_ascii(filename);
		delete ttecplot;
		tecplot=new clTecplot;
		sprintf(filename,"%s/test_unstructured_2zones_cc.dat",inputdir);
		tecplot->read(filename);
		ttecplot=new clTecplot;
		(*ttecplot)=(*tecplot);
		delete tecplot;
		sprintf(filename,"%s/output_test_copy_constructor_3_unstructured_2zones_cc.dat",outputdir);
		ttecplot->write_ascii(filename);
		delete ttecplot;
		tecplot=new clTecplot;
		sprintf(filename,"%s/test_unstructured_noconn.dat",inputdir);
		tecplot->read(filename);
		ttecplot=new clTecplot;
		(*ttecplot)=(*tecplot);
		delete tecplot;
		sprintf(filename,"%s/output_test_copy_constructor_3_unstructured_noconn.dat",outputdir);
		ttecplot->write_ascii(filename);
		delete ttecplot;
	}

	return 0;
};

