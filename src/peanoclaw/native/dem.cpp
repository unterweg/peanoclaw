#include"dem.h"

const float DEM::INVALID = -1e38f;

DEM::DEM(void) : m_data(NULL), m_boundary_data(NULL), m_bound_bottom(NULL), m_bound_top(NULL), m_bound_left(NULL), m_bound_right(NULL) {
	m_dimension[0]	= m_dimension[1]	= 0;
	m_LowerLeft[0]	= m_LowerLeft[1]	= m_LowerLeft[2] = INVALID;
	m_UpperRight[0]	= m_UpperRight[1]	= m_UpperRight[2]= INVALID;
	m_boundary_size = 0;
}

DEM::~DEM(void) {
	clear();
}

DEM::DEM(const DEM& other) : m_data(NULL), m_boundary_data(NULL), m_bound_bottom(NULL), m_bound_top(NULL), m_bound_left(NULL), m_bound_right(NULL) {
	m_dimension[0] = m_dimension[1] = 0;
	*this = other;
}

bool DEM::save(const std::string& name) const {
	FILE* fptr = NULL;
#ifdef _MSC_VER
	if (fopen_s(&fptr,name.c_str(),"wb")) return false;
#else 
	fptr = fopen(name.c_str(),"wb");
	if (fptr==NULL) return false;
#endif
	bool bRes = save(fptr);
	fclose(fptr);
	return bRes;
}

bool DEM::save(FILE* stream) const {
	if (stream==NULL) return false;
	int dims[3];
	dims[0] = int(m_dimension[0]);
	dims[1] = int(m_dimension[1]);
	dims[2] = int(m_boundary_size);
	if (fwrite(dims,sizeof(int)*3,1,stream)!=1)					return false;
	if (fwrite(m_LowerLeft,sizeof(double)*3,1,stream)!=1)		return false;
	if (fwrite(m_UpperRight,sizeof(double)*3,1,stream)!=1)		return false;
	
	const size_t chunk = (8<<20)/sizeof(float);	// 8MB chunks
	size_t size = nPixels();
	float* ptr = m_data;
	while (size>0) {
		size_t toWrite = std::min<size_t>(chunk,size);
		if (fwrite(ptr,sizeof(float)*toWrite,1,stream)!=1)		return false;
		size-=toWrite;
		ptr+=toWrite;
	}

	if (dims[2]>0) {
		float *ptr = m_boundary_data;
		size = nBoundaryPixels();
		while (size>0) {
			size_t toWrite = std::min<size_t>(chunk,size);
			if (fwrite(ptr,sizeof(float)*toWrite,1,stream)!=1)	return false;
			size-=toWrite;
			ptr+=toWrite;
		}
	}
	return true;
}

bool DEM::load(const std::string& name) {
	FILE* fptr = NULL;
#ifdef _MSC_VER
	if (fopen_s(&fptr,name.c_str(),"rb")) return false;
#else
	fptr = fopen(name.c_str(),"rb");
	if (fptr==NULL) return false;
#endif
	bool bRes = load(fptr);
	fclose(fptr);
	return bRes;
}

bool DEM::load(FILE* stream) {
	if (stream==NULL) return false;
	clear();
	int dims[3];
	if (fread(dims,sizeof(int)*3,1,stream)!=1)					return false;
	m_dimension[0] = dims[0];
	m_dimension[1] = dims[1];
	m_boundary_size= dims[2];
	if (fread(m_LowerLeft,sizeof(double)*3,1,stream)!=1)		return false;
	if (fread(m_UpperRight,sizeof(double)*3,1,stream)!=1)		return false;
	m_data = new float[nPixels()];
	if (!m_data) return false;
	const size_t chunk = (8<<20)/sizeof(float); // 8MB chunks
	size_t size = nPixels();
	float* ptr = m_data;
	while (size>0) {
		size_t toRead = std::min<size_t>(chunk,size);
		if (fread(ptr,sizeof(float)*toRead,1,stream)!=1)		return false;
		size-=toRead;
		ptr+=toRead;
	}
	if (m_boundary_size>0) {
		size_t size = nBoundaryPixels();
		m_boundary_data = new float[size];
		if (!m_boundary_data) return false;
		float* ptr = m_boundary_data;
		while (size>0) {
			size_t toRead = std::min<size_t>(chunk,size);
			if (fread(ptr,sizeof(float)*toRead,1,stream)!=1)	return false;
			size-=toRead;
			ptr+=toRead;
		}
		m_bound_bottom	=  m_boundary_data;
		m_bound_left	= &m_bound_bottom[m_boundary_size*(m_dimension[0]+2*m_boundary_size)];
		m_bound_right	= &m_bound_left[m_boundary_size*m_dimension[1]];
		m_bound_top		= &m_bound_right[m_boundary_size*m_dimension[1]];
	}
	return true;
}
