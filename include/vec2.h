
//
// a basic 2D vector class
//

class vec2 {
//-----------------------------------------------------------------------------
public:
  // constructors for none, one, and two given values
  vec2() : _x( 0 ), _y( 0 ) {};
  vec2( const double v ) : _x( v ), _y( v ) {};
  vec2( const double x, const double y ) : _x( x ), _y( y ) {};

  // coordinate access
  double& x() { return _x; };
  const double x() const { return _x; };
  double& y() { return _y; };
  const double y() const { return _y; };

  double& operator [] ( const int i ) {
    if ( i == 0 )
      return _x;
    return _y;
  };

  const double operator [] ( const int i ) const {
    if ( i == 0 )
      return _x;
    return _y;
  };

  // equality operator
  const bool operator == ( const vec2& p ) const {
    if ( _x == p.x() && _y == p.y() )
      return true;
    return false;
  };

  // negation operator
  const vec2 operator - () const {
    return vec2( -_x, -_y );
  };

  // addition
  const vec2 operator + ( const vec2& p ) const {
    return vec2( _x + p.x(), _y + p.y() );
  };

  vec2& operator += ( const vec2& p ) {
    _x += p.x(); _y += p.y();
    return *this;
  };

  // subtraction
  const vec2 operator - ( const vec2& p ) const {
    return vec2( _x - p.x(), _y - p.y() );
  };

  vec2& operator -= ( const vec2& p ) {
    _x -= p.x(); _y -= p.y();
    return *this;
  };

  // multiplication by a scalar
  const vec2 operator * ( const double w ) const {
    return vec2( _x * w, _y * w );
  };

  friend const vec2 operator * ( const double w, const vec2& p ) {
    return p * w;
  };

  vec2& operator *= ( const double w ) {
    _x *= w; _y *= w;
    return *this;
  };

  // division by a scalar
  const vec2 operator / ( const double w ) const {
    return vec2( _x / w, _y / w );
  };

  friend const vec2 operator / ( const double w, const vec2& p ) {
    return p / w;
  };

  vec2& operator /= ( const double w ) {
    _x /= w; _y /= w;
    return *this;
  };

  // dot product of two vectors
  const double operator * ( const vec2& p ) const {
    return ( _x * p.x() + _y * p.y() );
  };

  // Euclidean norm
  const double norm() const {
    return sqrt( _x * _x + _y * _y );
  };

  // print the vector
  friend std::ostream& operator << ( std::ostream& s, const vec2& p )  {
    s  << p.x() << " " << p.y();
    return s;
  };

//-----------------------------------------------------------------------------
private:
  // coordinates
  double _x, _y;
};
