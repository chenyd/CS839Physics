//#include "AnimatedTetrahedronMesh.h"
#include "FiniteElementMesh.h"

#include <Eigen/Dense>

#include <map>

template<class T>
struct LatticeMesh : public FiniteElementMesh<T>
{
    using Base = FiniteElementMesh<T>;

    // from AnimatedTetrahedonMesh
    using Base::m_meshElements;
    using Base::m_particleX;
    using Base::initializeUSD;
    using Base::initializeTopology;
    using Base::initializeParticles;
    using Vector3 = typename Base::Vector3;

    // from FiniteElementMesh
    using Base::initializeUndeformedConfiguration;
    using Base::m_stepEndTime;

    std::array<int, 3> m_cellSize; // dimensions in grid cells
    int m_a, m_b, m_c; // radius of sphere in grid cells
    T m_gridDX;
\
	const int m_pinchRadius;

    std::vector<std::array<int, 3>> m_activeCells; // Marks the "active" cells in the lattice
    std::map<std::array<int, 3>, int> m_activeNodes; // Maps the "active" nodes to their particle index

    std::vector<Vector3> m_particleUndeformedX;
    std::vector<int> m_leftHandleIndices;
    std::vector<int> m_rightHandleIndices;
    Vector3 m_leftHandleDisplacement;
    Vector3 m_rightHandleDisplacement;

    LatticeMesh()
        :Base(1.e2, 1., 4., .05), m_pinchRadius(1)
    {
        m_leftHandleDisplacement  = Vector3(.2, .1, 0.);
        m_rightHandleDisplacement = Vector3(-.2, -.1, 0.);
    }

    void initialize()
    {
        initializeUSD("Demo3D.usda");

        // Activate cells within a sphere of radius m_radius (in cells)

        for(int cell_i = 0; cell_i < m_cellSize[0]; cell_i++)
        for(int cell_j = 0; cell_j < m_cellSize[1]; cell_j++)
        for(int cell_k = 0; cell_k < m_cellSize[1]; cell_k++){

            // int r = (cell_i - m_cellSize[0]/2) * (cell_i - m_cellSize[0]/2) +
            //        (cell_j - m_cellSize[1]/2) * (cell_j - m_cellSize[1]/2) +
            //        (cell_k - m_cellSize[2]/2) * (cell_k - m_cellSize[2]/2);

            //if(r <= m_radius * m_radius)
            //    m_activeCells.push_back(std::array<int, 3>{cell_i, cell_j, cell_k});
        	if( cell_i - m_cellSize[0]/2 <= m_a && cell_i - m_cellSize[0]/2 >= -m_a)
        	if( cell_j - m_cellSize[0]/2 <= m_b && cell_j - m_cellSize[0]/2 >= -m_b)
        	if( cell_k - m_cellSize[0]/2 <= m_c && cell_k - m_cellSize[0]/2 >= -m_c)
        		m_activeCells.push_back(std::array<int, 3>{cell_i, cell_j, cell_k});
        }

        std::cout << "Created a model including " << m_activeCells.size() << " lattice cells" <<std::endl;

        // Create (uniquely numbered) particles at the node corners of active cells

        for(const auto& cell: m_activeCells){
            std::array<int, 3> node;
            for(node[0] = cell[0]; node[0] <= cell[0]+1; node[0]++)
            for(node[1] = cell[1]; node[1] <= cell[1]+1; node[1]++)
            for(node[2] = cell[2]; node[2] <= cell[2]+1; node[2]++){
                auto search = m_activeNodes.find(node);
                if(search == m_activeNodes.end()){ // Particle not yet created at this lattice node location -> make one
                    m_activeNodes.insert({node, m_particleX.size()});
                    m_particleX.emplace_back(m_gridDX * T(node[0]), m_gridDX * T(node[1]), m_gridDX * T(node[2]));
                }
            }
        }
        std::cout << "Model contains " << m_particleX.size() << " particles" << std::endl;

        // Make tetrahedra out of all active cells (6 tetrahedra per cell)

        for(const auto& cell: m_activeCells){
            int vertexIndices[2][2][2];
            for(int i = 0; i <= 1; i++)
            for(int j = 0; j <= 1; j++)
            for(int k = 0; k <= 1; k++){
                std::array<int, 3> node{cell[0] + i, cell[1] + j, cell[2] + k};
                auto search = m_activeNodes.find(node);
                if(search != m_activeNodes.end())
                    vertexIndices[i][j][k] = search->second;
                else
                    throw std::logic_error("particle at cell vertex not found");
            }

            m_meshElements.push_back(std::array<int, 4>{ vertexIndices[0][0][0], vertexIndices[1][0][0], vertexIndices[1][1][0], vertexIndices[1][1][1]});
            m_meshElements.push_back(std::array<int, 4>{ vertexIndices[0][0][0], vertexIndices[1][0][0], vertexIndices[1][1][1], vertexIndices[1][0][1]});
            m_meshElements.push_back(std::array<int, 4>{ vertexIndices[0][0][0], vertexIndices[1][0][1], vertexIndices[1][1][1], vertexIndices[0][0][1]});
            m_meshElements.push_back(std::array<int, 4>{ vertexIndices[0][0][0], vertexIndices[1][1][1], vertexIndices[0][1][1], vertexIndices[0][0][1]});
            m_meshElements.push_back(std::array<int, 4>{ vertexIndices[0][0][0], vertexIndices[1][1][1], vertexIndices[0][1][0], vertexIndices[0][1][1]});
            m_meshElements.push_back(std::array<int, 4>{ vertexIndices[0][0][0], vertexIndices[1][1][0], vertexIndices[0][1][0], vertexIndices[1][1][1]});
        }
        
        // Perform the USD-specific initialization of topology & particles
        // (this will also create a boundary *surface* to visualuze

        initializeTopology();
        initializeParticles();

        // Check particle indexing in mesh

        for(const auto& element: m_meshElements)
            for(const auto vertex: element)
                if(vertex < 0 || vertex >= m_particleX.size())
                    throw std::logic_error("mismatch between mesh vertex and particle array");
    
        // Initialize rest shape matrices and particle mass
        initializeUndeformedConfiguration();

        // Also record rest shape
        m_particleUndeformedX = m_particleX;

        // Identify particles on left and right handles

        for(const auto& cell: m_activeCells){
        	if(cell[0] - m_cellSize[0]/2 == m_a){
                std::array<int, 3> node{cell[0], cell[1], cell[2]};
                auto search = m_activeNodes.find(node);
                if(search != m_activeNodes.end())
                    m_leftHandleIndices.push_back(search->second);
                else
                    throw std::logic_error("particle at cell vertex not found");
            }

            if(cell[0] - m_cellSize[0]/2 == -m_a){
                std::array<int, 3> node{cell[0], cell[1], cell[2]};
                auto search = m_activeNodes.find(node);
                if(search != m_activeNodes.end())
                    m_rightHandleIndices.push_back(search->second);
                else
                    throw std::logic_error("particle at cell vertex not found");
            }
            
        }
    }

    void initializeDeformation()
    {
        // No need to apply any deformation; this example is driven by moving handles
    }

    void clearConstrainedParticles(std::vector<Vector3>& x) override
    { 
        for(const auto v: m_leftHandleIndices)
            x[v] = Vector3::Zero();
        for(const auto v: m_rightHandleIndices)
            x[v] = Vector3::Zero();       
    }

    void setBoundaryConditions() override
    { 
        T effectiveTime = std::min<T>(m_stepEndTime, 1.0);
        
        for(const auto v: m_leftHandleIndices){
            m_particleX[v] = m_particleUndeformedX[v] + effectiveTime * m_leftHandleDisplacement;
        }
        for(const auto v: m_rightHandleIndices){
            m_particleX[v] = m_particleUndeformedX[v] + effectiveTime * m_rightHandleDisplacement;
        }
    }
};

int main(int argc, char *argv[])
{
    LatticeMesh<float> simulationMesh;
    simulationMesh.m_cellSize = { 40, 40, 40 };
    simulationMesh.m_a = 8;
    simulationMesh.m_b = 4;
    simulationMesh.m_c = 4;
    simulationMesh.m_gridDX = 0.05;
    simulationMesh.m_nFrames = 10;
    simulationMesh.m_subSteps = 1;
    simulationMesh.m_frameDt = 0.1;

    // Initialize the simulation example
    simulationMesh.initialize();
    simulationMesh.initializeDeformation();
    // Output the initial shape of the mesh
    simulationMesh.writeFrame(0);

    // Perform the animation, output results at each frame
    for(int frame = 1; frame <= simulationMesh.m_nFrames; frame++){
        simulationMesh.simulateFrame(frame);
        simulationMesh.writeFrame(frame);
    }

    // Write the entire timeline to USD
    simulationMesh.writeUSD();

    return 0;
}

