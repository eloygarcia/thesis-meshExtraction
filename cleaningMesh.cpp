#include "cleaningMesh.h"

// public Constructor & Destructor;
cleaningMesh::cleaningMesh(void)
{
}

cleaningMesh::~cleaningMesh(void)
{
}

//
void cleaningMesh::removeDuplicatePoints(std::vector<double> &points_vector, std::vector<int> &cell_vector)
{
	std::cout << "removeDuplicatePoints" << std::endl;
	bool duplicate=false;
	std::vector<double> temp_log;
	
	double * temp_point_1 = new double[3]; // = {0.0,0.0,0.0};
	//	temp_point_1[0]=0.0;
	//	temp_point_1[1]=0.0;
	//	temp_point_1[2]=0.0;
	double * temp_point_2 = new double[3]; // = {0.0,0.0,0.0};
	//	temp_point_2[0]=0.0;
	//	temp_point_2[1]=0.0;
	//	temp_point_2[2]=0.0;

	for(int id = 0; id<(points_vector.size()/3); id++)
	{
		duplicate=false;
		
		temp_point_1[0]=points_vector[3*id];
		temp_point_1[1]=points_vector[3*id+1];
		temp_point_1[2]=points_vector[3*id+2];

		if(temp_point_1[0]==0 && temp_point_1[1]==0 && temp_point_1[2]==0)
		{
			duplicate = true;
		}

		for(int id2 = id+1; id2<(points_vector.size()/3); id2++)
		{
			if(id!=id2){
				temp_point_2[0]=points_vector[3*id2];
				temp_point_2[1]=points_vector[3*id2+1];
				temp_point_2[2]=points_vector[3*id2+2];

				if(temp_point_1[0]==temp_point_2[0] && temp_point_1[1]==temp_point_2[1] && temp_point_1[2]==temp_point_2[2])
				{
					duplicate = true;
					for(int z=0 ; z<cell_vector.size(); z++)
					{					
						if( cell_vector[z]==id ) cell_vector[z]=id2; 
					}
				}
			}
		}
		if(!duplicate)
		{
			temp_log.push_back( temp_point_1[0] );
			temp_log.push_back( temp_point_1[1] );
			temp_log.push_back( temp_point_1[2] );
		} 
	}
	 
	// =================== Return vector ======================
	points_vector.clear();
	points_vector = temp_log;

	// ===== Delete ====
	temp_log.~vector();
	delete [] temp_point_1;
	delete [] temp_point_2;
}

void cleaningMesh::removeDuplicatedCells( std::vector<int> &cell_vector)
{
	std::cout << "removeDuplicateCells" << std::endl;
	bool duplicate = false;
	std::vector<int> temp_cells;

	int * cell_1 = new int[3]; //= {0,0,0};
	//	cell_1[0]=0;
	//	cell_1[1]=0;
	//	cell_1[2]=0;
	int * cell_2 = new int[3]; // = {0,0,0};
	//	cell_2[0]=0;
	//	cell_2[1]=0;
	//	cell_2[2]=0;
	
	for(int i=0; i<(cell_vector.size()/3); i++)
	{
		duplicate = false; 

		cell_1[0]=cell_vector[3*i];
		cell_1[1]=cell_vector[3*i+1];
		cell_1[2]=cell_vector[3*i+2];
		
		if( cell_1[0]==cell_1[1] || cell_1[1]==cell_1[2] || cell_1[0]==cell_1[2])
		{
			duplicate = true;
		}

		for(int j=i+1; j<(cell_vector.size()/3); j++)
		{
			if(i!=j) 
			{
			cell_2[0]=cell_vector[3*j];
			cell_2[1]=cell_vector[3*j+1];
			cell_2[2]=cell_vector[3*j+2];

			if( cell_1[0]==cell_2[0] && cell_1[1]==cell_2[1] && cell_1[2]==cell_2[2]) duplicate=true;
			if( cell_1[0]==cell_2[0] && cell_1[1]==cell_2[2] && cell_1[2]==cell_2[1]) duplicate=true;
			if( cell_1[0]==cell_2[2] && cell_1[1]==cell_2[1] && cell_1[2]==cell_2[0]) duplicate=true;
			if( cell_1[0]==cell_2[1] && cell_1[1]==cell_2[0] && cell_1[2]==cell_2[2]) duplicate=true;
			if( cell_1[0]==cell_2[1] && cell_1[1]==cell_2[2] && cell_1[2]==cell_2[0]) duplicate=true;
			if( cell_1[0]==cell_2[2] && cell_1[1]==cell_2[0] && cell_1[2]==cell_2[1]) duplicate=true;
			}
		}

		if(!duplicate) // En este caso queda la última celda en lugar de la primera!
		{
			temp_cells.push_back(cell_1[0]);
			temp_cells.push_back(cell_1[1]);
			temp_cells.push_back(cell_1[2]);
		}
	}

	std::list<int> list;
	for(int m=0; m<temp_cells.size(); m++)
	{
		list.push_back( temp_cells[m] );
	}

	list.sort();
	list.unique();
//	list.sort();
	
	int count=0;
	//int max=0;
	for(std::list<int>::iterator id_list =list.begin(); id_list!=list.end(); id_list++)
	{
		for(int id_ce=0; id_ce<temp_cells.size(); id_ce++)
		{
			if( temp_cells[id_ce] == *id_list ) temp_cells[id_ce]=count; 
			//if( max<(*id_list) ) max= (*id_list);
		}
		count+=1;
	}
	
	// =================== return cell vector ====================
	cell_vector.clear();
	cell_vector = temp_cells;
	
	// ==== Delete ====
	temp_cells.~vector();

	delete [] cell_1;
	delete [] cell_2;
//
//	list.~list();
}

void cleaningMesh::searchCount( int k, int a, int b, std::vector<int> cell_vector,  int &ni , int &nd)
{
	ni=0; nd=0;
	int cell[3] = {0 ,0, 0};
	for(int m=0; m<cell_vector.size()/3; m++)
	{
		cell[0]=cell_vector[3*m];
		cell[1]=cell_vector[3*m+1];
		cell[2]=cell_vector[3*m+2];
		if(m!=k)
		{
			if( (cell[0]==a && cell[1]==b) || (cell[1]==a && cell[2]==b) || (cell[2]==a && cell[0]==b) ) nd+=1;
			if( (cell[1]==a && cell[0]==b) || (cell[2]==a && cell[1]==b) || (cell[0]==a && cell[2]==b) ) ni+=1;
		}
	}
}

void cleaningMesh::foundNonManifold(std::vector<int> cell_vector, std::list<int> &borderlist)
{
	std::cout << "foundNonManifold" << std::endl;
	int * cell = new int[3]; // ={0,0,0};
	//	cell[0] = 0;
	//	cell[1] = 0;
	//	cell[2] = 0;
	int ni_01 = 0; int nd_01 = 0;
	int ni_12 = 0; int nd_12 = 0;
	int ni_20 = 0; int nd_20 = 0;
	for(int n_c =0; n_c< cell_vector.size()/3; n_c++)
	{
		cell[0] =cell_vector[n_c*3];
		cell[1] =cell_vector[n_c*3+1];
		cell[2] =cell_vector[n_c*3+2];

		searchCount(n_c, cell[0], cell[1], cell_vector, ni_01, nd_01);
		searchCount(n_c, cell[1], cell[2], cell_vector, ni_12, nd_12);
		searchCount(n_c, cell[2], cell[0], cell_vector, ni_20, nd_20);

		if( (ni_01+nd_01==0) || ni_01+nd_01>1)
		{
			// borderlist.
			borderlist.push_back( cell[0] );
			borderlist.push_back( cell[1] );
		}
		if( (ni_12+nd_12==0) || ni_12+nd_12>1)
		{
			// borderlist.
			borderlist.push_back( cell[1] );
			borderlist.push_back( cell[2] );
		}
		if( (ni_20+nd_20==0) || ni_20+nd_20>1)
		{
			// borderlist.
			borderlist.push_back( cell[2] );
			borderlist.push_back( cell[0] );
		}
	}

	// ================ border list =====================
	borderlist.sort();
	borderlist.unique();

	// ==== Delete ====
	delete[] cell;
}

void cleaningMesh::performNewMesh( std::vector<double> &points_vector, std::vector<int> &cell_vector)
{
	std::list<int> clear_mesh;
	for(int v_c=0; v_c<cell_vector.size(); v_c++) clear_mesh.push_back( cell_vector[v_c] ); 

	clear_mesh.sort();
	clear_mesh.unique();

	std::vector<double> temp_points;
	std::vector<int> temp_cells;
	
	bool has_conexion;
	for(int id_point=0; id_point<(points_vector.size()/3); id_point++)
	{ 
		has_conexion=false;
		for(std::list<int>::iterator id_list =clear_mesh.begin(); id_list!=clear_mesh.end(); id_list++)
		{
			if( id_point== *id_list) has_conexion=true;
		}
		if( has_conexion) 
		{
			temp_points.push_back( points_vector[3*id_point]);
			temp_points.push_back( points_vector[3*id_point+1]);
			temp_points.push_back( points_vector[3*id_point+2]);
		}
	}

	int count=0;
	for(std::list<int>::iterator id_list =clear_mesh.begin(); id_list!=clear_mesh.end(); id_list++)
	{
		for(int id_ce=0; id_ce<cell_vector.size(); id_ce++)
		{
			if( cell_vector[id_ce] == *id_list ) cell_vector[id_ce]=count; 
			//if( max<(*id_list) ) max= (*id_list);
		}
		count+=1;
	}

	points_vector.clear();
	points_vector = temp_points;

	// ==== Delete =====
//	clear_mesh.~list();

	temp_points.~vector();
	temp_cells.~vector();
}

void cleaningMesh::clearMesh(std::vector<int> &cell_vector, std::list<int> borderlist, bool &repeat)
{
	std::cout << "clearMesh" << std::endl;
	bool repeat_1 = false; bool repeat_2 = false; bool repeat_3 = false;
	std::vector<int> temp_cell;
	int cell[3] = {0, 0, 0};
	for(int v_c=0; v_c<cell_vector.size()/3; v_c++)
	{
		repeat = false;
		bool is_1 = false; bool is_2 = false; bool is_3 = false;
		cell[0]=cell_vector[3*v_c];
		cell[1]=cell_vector[3*v_c+1];
		cell[2]=cell_vector[3*v_c+2];

		for(std::list<int>::iterator v_list=borderlist.begin(); v_list!=borderlist.end(); v_list++)
		{
			if( cell[0] == *v_list ) is_1 = true;
			if( cell[1] == *v_list ) is_2 = true;
			if( cell[2] == *v_list ) is_3 = true;
		}

		if( !(is_1 && is_2 && is_3))
		{
			temp_cell.push_back( cell[0] );
			temp_cell.push_back( cell[1] );
			temp_cell.push_back( cell[2] );
		}
	}

	repeat = false;
	// ============================= Return Cell_vector =================================
	cell_vector.clear();
	cell_vector = temp_cell; 

	// ==== Delete ====
	temp_cell.~vector();
}