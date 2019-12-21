#include "convertfvtotecplot.h"

#ifdef DEBUG
template <> int clVector_namespace::clVector<int>::_num_vector=0;
template <> int clVector_namespace::clVector<float>::_num_vector=0;
template <> int clVector_namespace::clVector<double>::_num_vector=0;
template <> int clVector_namespace::clMatrix<int>::_num_matrix=0;
template <> int clVector_namespace::clMatrix<double>::_num_matrix=0;
template <> int clVector_namespace::clMatrix3D<double>::_num_matrix=0;
int clString_namespace::clString::_num_string=0;
int clString_namespace::clString_vector::_num_string_vector=0;
#endif

int main(int argc, char* argv[])
{
	using clUtils_namespace::fatal_err;
	using clUtils_namespace::read_int;

	int num_file=0,filetype=0;
	const int NAME_LENGTH=512;
	char gridname[NAME_LENGTH],filename[NAME_LENGTH];
	if (argc<=1) {
		num_file=2;
		fprintf(stdout,"Enter file type (1: single file, 2: separate grid/solution files): ");
		fscanf(stdin,"%d",&num_file);
		if (num_file==1) {
			fprintf(stdout,"Enter FieldView filename: ");
			fscanf(stdin,"%s",filename);
		} else if (num_file==2) {
			fprintf(stdout,"Enter FieldView grid filename: ");
			fscanf(stdin,"%s",gridname);
			fprintf(stdout,"Enter FieldView solution filename: ");
			fscanf(stdin,"%s",filename);
		} else {
			fatal_err("Invalid file type");
		}
		fprintf(stdout,"Output binary (0) or ASCII (1) Tecplot files: ");
		fscanf(stdin,"%d",&filetype);
	} else if (argc==2) {
		num_file=1;
		strcpy(filename,argv[1]);
		fprintf(stdout,"Reading FieldView filename: %s\n",filename);
		filetype=0;
	} else if (argc==3) {
		num_file=2;
		strcpy(filename,argv[1]);
		fprintf(stdout,"Reading FieldView filename: %s\n",gridname);
		strcpy(filename,argv[2]);
		fprintf(stdout,"Reading FieldView filename: %s\n",filename);
		filetype=0;
	} else {
		fatal_err("Invalid number of command line arguments in main");
	}

	/* Read a binary FieldView file(s) */
	clFieldview *fieldviewgrid=NULL;
	if (num_file==2) {
		fieldviewgrid=new clFieldview;
		fieldviewgrid->read(gridname);
	}
	clFieldview *fieldview=new clFieldview;
	fieldview->read(filename);

	/* Convert to Tecplot and print */
	clTecplot *tecplot=new clTecplot;
	if (num_file==1) {
		fieldview->convert_to_tecplot_surface(tecplot);
	} else {
		fieldviewgrid->convert_to_tecplot_surface(tecplot);
		fieldview->add_result_to_tecplot_surface(tecplot);
	}

	const int num_zone=tecplot->get_num_zone();
	clVector_namespace::clVector<int> writezone(num_zone,1);
	FILE *fileptr=fopen("select_zones.dat","r");
	if (fileptr) {
		int n_zone,n_write;
		const int num_write=read_int(fileptr);
		if (num_write<0 || num_write>=num_zone) fatal_err("Illegal value for num_zone");
		writezone=0;
		for (n_write=0; n_write<num_write; n_write++) {
			n_zone=read_int(fileptr)-1;
			if (n_zone<0 || n_zone>=num_zone)
				clUtils_namespace::warning("Illegal zone value in select_zone.dat");
			writezone(n_zone)=1;
		}
		fclose(fileptr);
		int shareflag=0;
		for (n_zone=0; n_zone<num_zone; n_zone++)
			if (!writezone(n_zone) && tecplot->get_zonesharing(n_zone)) shareflag++;
		if (shareflag) tecplot->eliminate_varsharing();
	}

	const int num_var=tecplot->get_num_var();
	clVector_namespace::clVector<int> writevar(num_var,1);
	fileptr=fopen("select_variables.dat","r");
	if (fileptr) {
		int n_var,n_write,varflag;
		const int num_write=read_int(fileptr);
		if (num_write<0 || num_write>=num_var) fatal_err("Illegal value for num_write");
		clString_namespace::clString varname,vvarname;
		writevar(n_var)=0;
		for (n_write=0; n_write<num_write; n_write++) {
			varname.read_ascii(fileptr);
			++varname;
			varflag=0;
			for (n_var=0; n_var<num_var; n_var++) {
				vvarname=tecplot->get_varname(n_var);
				++vvarname;
				if (varname==vvarname) {
					writevar(n_var)=1;
					varflag++;
					break;
				}
			}
			if (!varflag) fatal_err("Unable to match variable name");
		}
		fclose(fileptr);
	}

	/* Print Tecplot file */
	if (filetype) {
		fprintf(stdout,"Outputing ASCII Tecplot file\n");
		clUtils_namespace::fileextension(filename,"dat",NAME_LENGTH);
		tecplot->write_ascii(writezone,filename,writevar);
	} else {
		fprintf(stdout,"Outputing binary Tecplot file\n");
		clUtils_namespace::fileextension(filename,"plt",NAME_LENGTH);
		tecplot->write(filename,writevar);
	}

	if (num_file==2) delete fieldviewgrid;
	delete fieldview;
	delete tecplot;

	return 0;
}

