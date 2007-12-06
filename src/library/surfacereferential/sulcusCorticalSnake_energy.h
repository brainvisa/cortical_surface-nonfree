#ifndef AIMS_CORTICAL_SNAKE_ENERGY_H
#define AIMS_CORTICAL_SNAKE_ENERGY_H


namespace aims
{

	class SulcusCorticalSnake_energy
	{
	
		public:
			
			//Global variables
			//Update when class is instanciated
			
			//vector list_points = liste des indices des points Pi decrivant le snake
			std::vector<uint> list_points;
			
			//index du sommet changeant (DANS LA LISTE list_points!)
			uint index_courant;
			
			uint n1, n2;
			
			float alpha1, alpha2, alpha3;
			
			AimsSurfaceTriangle mesh;
			
			AimsSurface<3, Void> mesh_base;
			
			uint size;
			
			float h_min, h_max;
			
			TimeTexture<float> curv;
			
			float max;			
			
			TimeTexture<float> tex_distance;
			
			std::vector<Point3df> vert;
			
			std::vector< AimsVector<uint,3> > poly;
			
			uint size_bucket;
			
			uint size_vector;
			
			
			//Constructor
			
			SulcusCorticalSnake_energy( std::vector<uint> _list_points, uint _index_courant, uint _n1, uint _n2, float _alpha1, float _alpha2, float _alpha3, AimsSurfaceTriangle _mesh, uint _size, float _h_min, float _h_max, TimeTexture<float> _curv, float _max, TimeTexture<float> _tex_distance ) : list_points(_list_points), index_courant(_index_courant), n1(_n1), n2(_n2), alpha1(_alpha1), alpha2(_alpha2), alpha3(_alpha3), mesh(_mesh), size(_size), h_min(_h_min), h_max(_h_max), curv(_curv), max(_max), tex_distance(_tex_distance)
			{
// 				std::cout<<"Constructeur SulcusCorticalSnake_energy"<<std::endl;
				size_vector=list_points.size();
			//Mesh version AimsSurface
				mesh_base=mesh[0];
				vert = mesh_base.vertex();
				poly=mesh.polygon();
				
			}
			
			~SulcusCorticalSnake_energy()
			{
// 				std::cout<<"Destructeur SulcusCorticalSnake_energy"<<std::endl;
				vert.clear();
				poly.clear();
				list_points.clear();
				
				tex_distance.erase();
				curv.erase();
			}
			
			//Functions
			
			float total_energy();
// 			float distance_energy();
// 			float distance_Pi_to_bucket(uint index_vector);
// 			float weighted_square_distance(uint i, Point3d j, float weight);
// 			float square_distance(uint i, Point3d j);
			float curvature_energy();
			float elastic_energy();
			std::vector<float> geodesic_distance(uint origin, uint i, uint j);
			std::vector<float> MeshDistance_adapt( const Texture<float> & inittex, uint uind1, uint ind2 );
	};
	
}



#endif










