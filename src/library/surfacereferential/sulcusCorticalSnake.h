#ifndef AIMS_SULCUS_CORTICAL_SNAKE_H
#define AIMS_SNAKE_H


#include <aims/mesh/texture.h>
#include <aims/io/reader.h>
#include <aims/io/writer.h>
#include <aims/mesh/surfaceOperation.h>
#include <cortical_surface/surfacereferential/sulcusCorticalSnake_energy.h>

namespace aims
{

	class SulcusCorticalSnake
	{
	
		public:
			
			//Global variables
			//Update when class is instanciated
			TimeTexture<float> constraint_load;
			
			float value;
			float alpha1, alpha2, alpha3;
			AimsSurfaceTriangle mesh;
			
			uint n1,n2, size;
			uint size_vector;
			int add_pts;
			int cpt_multi;
			TimeTexture<float> curv;
			
			float h_min, h_max;
			
			float max;
			
			TimeTexture<float> constraint;
			TimeTexture<float> result_total;
			TimeTexture<float> tex_distance;
			
		//Liste de voisinage
			std::vector<std::set<uint> > neigh;
			
			//vector list_points = liste des indices des points Pi decrivant le snake
			std::vector<uint> list_points;
			
			//vector list_points = liste des indices, pour comparer et mettre a jour cpt_points
			std::vector<uint> avant_list;
			
			//vector cpt_points = compteur des points Pi decrivant le snake
			std::vector<uint> cpt_points;
			
			
			//vector ordre_points = position du point dans le snake
			std::vector<uint> ordre_points;
			
			//vector new_vector = liste des indices des points Pi decrivant le snake pour comparaison pour stop_global
			std::vector<uint> new_vector;
			
			//vector new_vector_res = liste des indices des points Pi decrivant le snake pour une resolution donnee
			std::vector<uint> new_vector_res;
			
			std::vector<Point3df> vert;
			AimsSurface<3, Void> mesh_base;
			std::vector< AimsVector<uint,3> > poly;
			
			int cpt_resolution;
			
// 			TimeTexture<float> res;
			
			//Constructor
			
			SulcusCorticalSnake( TimeTexture<float> _constraint, float _value, float _alpha1, float _alpha2, float _alpha3, AimsSurfaceTriangle _mesh) : constraint_load(_constraint), value(_value), alpha1(_alpha1), alpha2(_alpha2), alpha3(_alpha3), mesh(_mesh)
			{
				//Computing neighbourhood
// 				std::cout<<"Constructor!"<<std::endl;
				//			std::cout<<"Constructeur SulcusCorticalSnake"<<std::endl;
				neigh = SurfaceManip::surfaceNeighbours( mesh );
				mesh_base=mesh[0];
				vert = mesh_base.vertex();
				poly = mesh_base.polygon();
				
				size=constraint_load[0].nItem();
				
// 				std::cout<<"value="<<value<<" - alpha1="<<alpha1<<std::endl;
				
				for(uint i=0;i<size;i++)
				{
					if(fabs(constraint_load[0].item(i)-_value) < 0.01)
						constraint[0].push_back(_value);
					else
						constraint[0].push_back(0.0);
				}
// 				std::cout<<"size copie"<<constraint[0].nItem()<<std::endl;
				for( uint i=0; i<constraint[0].nItem(); i++)
				{
					result_total[0].push_back(0);
				}
				cpt_multi=0;
				cpt_resolution=0;
				
				/*
				Writer<Texture1d> wr1fd("cartea.tex");
				wr1fd.write(constraint);*/
			}
			
			~SulcusCorticalSnake()
			{
				//			std::cout<<"Destructeur SulcusCorticalSnake"<<std::endl;
				neigh.clear();
				vert.clear();
				poly.clear();
				list_points.clear();
				avant_list.clear();
				cpt_points.clear();
				ordre_points.clear();
				new_vector.clear();
				new_vector_res.clear();
				
				constraint_load.erase();
				constraint.erase();
				result_total.erase();
				tex_distance.erase();
			}
			
			//Functions
			void createDistanceTexture();
			TimeTexture<float> compute_snake();
// 			void compute_snake_at_1_resolution(TimeTexture<float> & result_multi);
			void compute_snake_at_1_resolution();
			int stop_total();
			int stop_condition_1_resolution();
			int define_extremities();
			int is_it_in_the_vector(uint vertex);
			int are_they_neighbours(uint one, uint two);
			uint define_new_middle_point(uint i, uint j);
			void refine_vector();
			void compute_curv();
			float compute_energy(uint index_courant);
			void treat_list_point(uint index, std::map< int, std::map< int, int > > & count);
			void process_list( std::map< int, std::map< int, int > > & count );
			float MeshDistance_adapt( const Texture<float> & inittex, bool allowUnreached, uint ind );
			Texture<float> MeshDistance_adapt_tex( const Texture<float> & inittex, bool allowUnreached, uint ind );
			std::map< uint, float> MeshDistance_adapt_local( const Texture<float> & inittex, uint ind1, uint ind2 );
	};
	
}


#endif
