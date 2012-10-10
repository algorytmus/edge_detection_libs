
#include <math.h>
#include <cmath>
#include <conio.h>
#include <stdio.h>
#include <limits>

using namespace std;

/*
Pawel Pietrzak, p.d.pietrzak@gmail.com
Komputerowa analiza i przetwarzanie obrazów cyfrowych - projekt
*/


template <int Size,class IM>
class Square
{
public:
	typename typedef IM::ValueType ValueType;
	static const int SquareSize=Size;
	typename typedef IM::CoordinateType CoordinateType;
	virtual inline ValueType& operator()(CoordinateType x,CoordinateType y)=0;
};

template <typename DirectionType>
class Direction;
/*Pochodna dwuwymiarowa*/
template <typename ValueType>
struct Derivative
{
	ValueType Gy;
	ValueType Gx;
};
/*Klasa zakresu dla zmiennych rozmiarów okien rozwiajana statycznie*/
template <int Size>
class Range
{
public:
	static const int left=((Size%2==0?((-(Size/2))+1):(-(Size/2))));
	static const int right=Size/2;
	static const int offset=-left;
};

template <int Size,typename ValueType,typename CoordinateType,ValueType (*fun)(CoordinateType,CoordinateType)>
class Mask
{

public:
	typename typedef ValueType MaskType;
	typedef ValueType (*FunctionType)(CoordinateType,CoordinateType);
	static const FunctionType Function;
public:
	static inline ValueType mask(CoordinateType x,CoordinateType y)
	{
		return fun(x,y);
	}
	static const ValueType tab[Size][Size];
private:
	static const bool dummy_field;
};
template <int Size,typename ValueType,typename CoordinateType,ValueType (*fun)(CoordinateType,CoordinateType)>
const bool Mask<Size,ValueType,CoordinateType,fun>::Function=fun;

template <class Mask,class Matrix, class SumAccumulatorType>
void filter(Matrix& matrix)
{	
	/*
		Dyskretna wersja G(x,y,sigma)
		sigma=1.4
	*/
	typename typedef Matrix::ValueType ValueType;
	const int Size= Matrix::SquareSize;
	SumAccumulatorType newValue=SumAccumulatorType();
	for(int i=Range<Size>::left;i<=Range<Size>::right;i++)
			for(int j=Range<Size>::left;j<=Range<Size>::right;j++)
				newValue+=(matrix(j,i))*(Mask::mask(j,i));
	matrix(0,0)=newValue/159;
}

template <class Mask,class Matrix, class SumAccumulatorType>
void filter_log(Matrix& matrix,Matrix& outputMatrix)
{	
	typename typedef Matrix::ValueType MatrixType;
	typename typedef Mask::MaskType MaskType;
	const int Size= Matrix::SquareSize;
	SumAccumulatorType newValue=SumAccumulatorType();
	for(int i=Range<Size>::left;i<=Range<Size>::right;i++)
			for(int j=Range<Size>::left;j<=Range<Size>::right;j++)
			{
				newValue+=((MatrixType)matrix(j,i))*(MaskType)Mask::mask(j,i);
			}
	outputMatrix(0,0)=(MatrixType)newValue;
}



template <int Size,class IM>
class Square;
/*Abstrakcyjna klasa pozwalaj¹ca na iterowanie okien o zadanym statycznie rozmiarze*/
template <typename ValueType,typename CoordinateType>
class ImageMapper
{

public:
	typename typedef ValueType ValueType;
	typename typedef CoordinateType CoordinateType;
	typename typedef Derivative<ValueType> DerivativeType;
	typename typedef ImageMapper ThisType;
	virtual inline ValueType& getPixel(CoordinateType x,CoordinateType y)=0;
	virtual inline CoordinateType getLeftTopX()=0;
	virtual inline CoordinateType getLeftTopY()=0;
	virtual inline CoordinateType getRightBottomX()=0;
	virtual inline CoordinateType getRightBottomY()=0;
	template <int Size>
	class Iterator
      {
                    public:
						typename typedef ValueType ValueType;
					typename typedef Square<Size,ThisType> SquareType;
						class IterSquare: public Square<Size,ThisType>
						{
							public:
							virtual ValueType& operator()(CoordinateType x,CoordinateType y)
							{
								return t->getPixel((*_x)+x,(*_y)+y);
							}
							virtual ValueType& get(CoordinateType x,CoordinateType y)
							{
								return t->getPixel((*_x)+x,(*_y)+y);
							}
							template <typename DirectionValue>
							inline bool canTravel(DirectionValue direction)
							{
								return it->canTravel<DirectionValue>(direction);
							}
							template <typename DirectionValue>
							inline void travel(DirectionValue direction)
							{
								it->travel(direction);
							}
							CoordinateType* _x;
							CoordinateType* _y;
							ThisType* t;
							Iterator* it;
						};

						Iterator(ThisType* t1):t(t1)
						{
							x = t->getLeftTopX()+Range<Size>::offset;
							y = t->getLeftTopY()+Range<Size>::offset;
							endReached=false;
							square._x=&x;
							square._y=&y;
							square.t=t;
							square.it=this;
						}
						Iterator()
						{
							endReached=true;
						}
					template <typename DirectionValue>
					inline bool canTravel(const DirectionValue& direction)
					{
						CoordinateType xCor=x+Direction<DirectionValue>::neighborX<CoordinateType>(direction);
						if(xCor-Range<Size>::offset<t->getLeftTopX() || xCor+Range<Size>::right>t->getRightBottomX())
							return false;
						CoordinateType yCor=y+Direction<DirectionValue>::neighborY<CoordinateType>(direction);
						if(yCor-Range<Size>::offset<t->getLeftTopY() || yCor+Range<Size>::right>t->getRightBottomY())
							return false;
						return true;
					}
                   
					template <typename DirectionValue>
					inline void travel(DirectionValue direction)
					{
						y+=Direction<DirectionValue>::neighborY<CoordinateType>(direction);
						x+=Direction<DirectionValue>::neighborX<CoordinateType>(direction);
					}

					IterSquare& operator*()
                   {
                         return square;
                   }
                   Iterator& operator++(int)
                   {
						if((y+1+Range<Size>::right)>t->getRightBottomY())
					   {
                         endReached=true;
						return *this;
					   }
					   if((x+1+Range<Size>::right)>t->getRightBottomX())
					   {
						   x = t->getLeftTopX()+Range<Size>::offset;
						   y++;
					   }
					   else
					   {
						   x++;
					   }
					   
                         return *this;
                   }
                   bool operator==( const Iterator & iterator ) const
                   {
                         return endReached==iterator.endReached;
                   }
                   bool operator!=( const Iterator & iterator ) const
                   {
                         return endReached!=iterator.endReached;
                   }
           // private:
			
			
			CoordinateType x;
			CoordinateType y;
			bool endReached;
            IterSquare square; 
			ThisType* t;
      };
		template <int Size>
		Iterator<Size> begin()
		{
			return Iterator<Size>(this);
		}
		template <int Size>
		Iterator<Size> end()
		{
			return Iterator<Size>();
		}
		
	inline CoordinateType size()
	{
		return (getRightBottomX()-getLeftTopX()+1)*(getRightBottomY()-getLeftTopY()+1);
	}


	void init()
	{
		for(int y=getLeftTopY();y<=getRightBottomY();y++)
			for(int x=getLeftTopX();x<=getRightBottomX();x++)
				getPixel(x,y)=ValueType();
	}
	protected:
	
	inline virtual void* mem()=0;
};



template <class Square,typename ResultType>
class DerivativePolicy
{
public:
	template <int Size>
	class DerivativeConcretePolicy
	{
		public:
		static inline Derivative<ResultType> getDerivativeInside(Square& square);
	};
	template <>
	class DerivativeConcretePolicy<2>
	{
		public:
		static inline Derivative<ResultType> getDerivativeInside(Square& square)
		{
			Derivative<ResultType> der;
			der.Gx=(square(1,0)-square(0,0))+(square(1,1)-square(0,1));
			der.Gy=(square(1,0)-square(1,1))+(square(0,0)-square(0,1));
			return der;
		}
	};
	template <>
	class DerivativeConcretePolicy<3>//drugiego rzêdu
	{
		public:
		static inline Derivative<ResultType> getDerivativeInside(Square& square)
		{
			Derivative<ResultType> der;
			der.Gx=square(-1,0)-(2*square(0,0))+square(1,0);
			der.Gy=square(0,1)-(2*square(0,0))+square(0,-1);
			return der;
		}
	};
	static inline Derivative<ResultType> getDerivative(Square& square)
	{
		return DerivativeConcretePolicy<Square::SquareSize>::getDerivativeInside(square);
	}
	
};





template <typename ValueType,typename CoordinateType>
class AllocatedValueMap
	: public ImageMapper<ValueType,CoordinateType>
	{
	public:
	AllocatedValueMap(CoordinateType _width,CoordinateType _height)
	{
		x_start=0;
		y_start=0;
		x_end = _width-1;
		y_end = _height-1;
		iSrc=new ValueType[(x_end-x_start+1)*(y_end-y_start+1)];
	}
	template <class IM>
	AllocatedValueMap(IM& im)
	{
		copySizeAndAllocate(im);
	}
	template <class IM>
	AllocatedValueMap(IM& im,bool copy)
	{
		copySizeAndAllocate(im);
		typename typedef IM::ValueType IMValueType;
		if(copy)
			mcopy<IMValueType>(im.mem());
	}
	
	virtual ValueType& getPixel(CoordinateType x,CoordinateType y)
	{
		return iSrc[y*(x_end-x_start+1)+x];
	}
	
	virtual inline CoordinateType getLeftTopX()
	{
		return x_start;
	}
	virtual inline CoordinateType getLeftTopY()
	{
		return y_start;
	}
	virtual inline CoordinateType getRightBottomX()
	{
		return x_end;
	}
	virtual inline CoordinateType getRightBottomY()
	{
		return y_end;
	}
	virtual ~AllocatedValueMap()
	{
		delete[] iSrc;
	}
	protected:
		inline void* mem()
		{
			return (void*)iSrc;
		}
	private:
	template <typename CopyType>
	inline void mcopy(void* m)
	{
		iSrc=new ValueType[size()];
		CopyType* ctM=(CopyType*)m;
		for(int i=0;i<size();i++)
			iSrc[i]=(ValueType)ctM[i];
	}
	template <>
	inline void mcopy<ValueType>(void* m)
	{
		memcpy(iSrc,m,size()*sizeof(ValueType));
	}
	template <class IM>
	void copySizeAndAllocate(IM& im)
	{
		x_start=im.getLeftTopX();
		y_start=im.getLeftTopY();
		x_end=im.getRightBottomX();
		y_end=im.getRightBottomY();
		iSrc=new ValueType[(x_end-x_start+1)*(y_end-y_start+1)];
	}
	CoordinateType x_start;
	CoordinateType y_start;
	CoordinateType x_end;
	CoordinateType y_end;
	ValueType *iSrc;
};




template <typename DirectionType>
class Direction
{
public:
	static const DirectionType none=0;
	static const DirectionType dir1=1;
	static const DirectionType dir2=2;
	static const DirectionType dir3=3;
	static const DirectionType dir4=4;
	template <class Square>
	static inline typename Square::ValueType& neighbor(Square& s,const DirectionType& dir)
	{
		return s(neighborX<Square::CoordinateType>(dir),neighborY<Square::CoordinateType>(dir));
	}
	template <typename CoordinateType>
	static inline CoordinateType neighborX(const DirectionType& dir)
	{
		switch(dir)
		{
			case dir1: 
				return 0;
			case dir2: 
				return 1;
			case dir3: 
				return 1;
			case dir4: 
				return 1;
			case -dir1: 
				return 0;
			case -dir2: 
				return -1;
			case -dir3: 
				return -1;
			case -dir4: 
				return -1;
		}
	}
	template <typename CoordinateType>
	static inline CoordinateType neighborY(const DirectionType& dir)
	{
		switch(dir)
		{
			case dir1: 
				return -1;
			case dir2: 
				return -1;
			case dir3: 
				return 0;
			case dir4: 
				return 1;
			case -dir1: 
				return 1;
			case -dir2: 
				return 1;
			case -dir3: 
				return 0;
			case -dir4: 
				return -1;
		}
	}
};

//NonMaximumSuppressionAlgorithms
class NonMaximumSuppressionAlgorithms
{
public:
	/*Kwadrat powinien byæ 3x3*/
	template <class Square>
	static inline typename Square::ValueType centerPixelBuildsEdge(Square& s)
	{
		typename typedef Square::ValueType ValueType;
		if(Direction<Square::ValueType>::dir1==s(0,-1) && Direction<Square::ValueType>::dir1==s(0,1))return Direction<Square::ValueType>::dir1;
		if(Direction<Square::ValueType>::dir3==s(-1,0) && Direction<Square::ValueType>::dir3==s(1,0))return Direction<Square::ValueType>::dir3;
		if(Direction<Square::ValueType>::dir4==s(-1,-1)&& Direction<Square::ValueType>::dir4==s(1,1))return Direction<Square::ValueType>::dir4;
		if(Direction<Square::ValueType>::dir2==s(-1,1) && Direction<Square::ValueType>::dir2==s(1,-1))return Direction<Square::ValueType>::dir2;
		return Direction<Square::ValueType>::none;
	}
	template <class AngleType,typename DirectionType>
	static inline typename DirectionType direction(AngleType val)
	{
			if ( ( (val <= 0.3927) && (val >= -0.3927) ) )
				return Direction<DirectionType>::dir1;//ok
			if((val >= -1.177) && (val <= -0.3927))
				return Direction<DirectionType>::dir2;//ok
			if( (val <= 1.177) && (val >= 0.3927) )
				return Direction<DirectionType>::dir4;//ok
			if ( ( (val >= 1.177) || (val <=-1.177) ))
				return Direction<DirectionType>::dir3;//ok
			return Direction<DirectionType>::none;
	}
	/*Kwadrat powinien byæ 3x3; Square s - kwadrat modu³ów; DirectionType inEdgeOfDirection - kierunek krawêdzi do jakiej nale¿y*/
	
	
	template <typename DirectionType,class Square>
	
	static inline bool classifiesToRemoval(Square& s,DirectionType inEdgeOfDirection)
	{
		typename typedef Square::ValueType ValueType;
		//if(inEdgeOfDirection!=Direction<DirectionType>::none)
		{
			ValueType midVal=s(0,0);
			switch(inEdgeOfDirection)
			{
			case Direction<DirectionType>::dir1:
					if(s(-1,0)>midVal || s(1,0)>midVal)
						return true;
					break;
			case Direction<DirectionType>::dir2:
					if(s(-1,-1)>midVal || s(1,1)>midVal)
						return true;
					break;
			case Direction<DirectionType>::dir3:
					if(s(0,-1)>midVal || s(0,1)>midVal)
						return true;
					break;
			case Direction<DirectionType>::dir4:
					if(s(-1,1)>midVal || s(1,-1)>midVal)
						return true;
					break;
			}
			return false;
		}
	}
};

class ZeroCrossing
{
public:
	template <class Square>
	static inline bool crossesZero2x2(Square& s,typename Square::ValueType slopeTresholdA,typename Square::ValueType slopeTresholdB)
	{
		typename typedef Square::ValueType ValueType;
		bool qual=false;
		if(s(0,0)>0 & s(1,0)<0) qual=true;
		if(s(0,0)<0 & s(1,0)>0) qual=true;
		if(s(0,0)>0 & s(0,1)<0) qual=true;
		if(s(0,0)<0 & s(0,1)>0) qual=true;
		if(s(0,0)<0 & s(1,1)>0) qual=true;
		if(s(0,0)>0 & s(1,1)<0) qual=true;

		if(qual)
			if(s(0,0)<slopeTresholdA || s(0,0)>slopeTresholdB)
				return true;
		return false;
	}
	template <class Square>
	static inline bool crossesZero3x3(Square& s,typename Square::ValueType slopeTresholdA,typename Square::ValueType slopeTresholdB)
	{
		typename typedef Square::ValueType ValueType;
		bool qual=false;
		if(s(-1,-1)>0 & s(1,1)<0) qual=true;
		if(s(-1,-1)<0 & s(1,1)>0) qual=true;

		if(s(0,-1)<0 & s(0,1)>0) qual=true;
		if(s(0,-1)>0 & s(0,1)<0) qual=true;

		if(s(1,-1)>0 & s(-1,1)<0) qual=true;
		if(s(1,-1)<0 & s(-1,1)>0) qual=true;

		if(s(-1,0)>0 & s(1,0)<0) qual=true;
		if(s(-1,0)<0 & s(1,0)>0) qual=true;

		if(qual)
			if(s(0,0)<slopeTresholdA || s(0,0)>slopeTresholdB)
				return true;
		return false;
	}

};

//end NonMaximumSuppressionAlgorithms

/////////HLvalues
template <typename TType>
class HLvalues
{
public:
	TType high;
	TType low;
};
////end HLvalues

/////HysteresisThresholding
template <typename DirectionType>
class HysteresisThresholding
{
private:
	template <class RemovalSquare,class ModuleSquare,DirectionType dir>
static void reasoninInDirection(ModuleSquare& modStr,RemovalSquare& rem,RemovalSquare & possibleToBePossitive,const HLvalues<typename ModuleSquare::ValueType>& vals)
	{
		if(Direction<DirectionType>::neighbor(rem,dir)
			&& !Direction<DirectionType>::neighbor(possibleToBePossitive,dir))
			if(Direction<DirectionType>::neighbor(modStr,dir)>vals.low)
			{
				Direction<DirectionType>::neighbor<RemovalSquare>(rem,dir)=false;
				if(rem.canTravel(dir))
				{
					rem.travel(dir);
					modStr.travel(dir);
					possibleToBePossitive.travel(dir);
					reasonPosotiveNeighbors(modStr,rem,possibleToBePossitive,vals);
					rem.travel(-dir);
					modStr.travel(-dir);
					possibleToBePossitive.travel(-dir);
				}
			}
	}
public:
template <class RemovalSquare,class ModuleSquare>
static void reasonPosotiveNeighbors(ModuleSquare& modStr,RemovalSquare& rem,RemovalSquare & possibleToBePossitive,const HLvalues<typename ModuleSquare::ValueType>& vals)
	{
		reasoninInDirection<RemovalSquare,ModuleSquare,Direction<DirectionType>::dir1>(modStr,rem,possibleToBePossitive,vals);
		reasoninInDirection<RemovalSquare,ModuleSquare,Direction<DirectionType>::dir2>(modStr,rem,possibleToBePossitive,vals);
		reasoninInDirection<RemovalSquare,ModuleSquare,Direction<DirectionType>::dir3>(modStr,rem,possibleToBePossitive,vals);
		reasoninInDirection<RemovalSquare,ModuleSquare,Direction<DirectionType>::dir4>(modStr,rem,possibleToBePossitive,vals);
		reasoninInDirection<RemovalSquare,ModuleSquare,-Direction<DirectionType>::dir1>(modStr,rem,possibleToBePossitive,vals);
		reasoninInDirection<RemovalSquare,ModuleSquare,-Direction<DirectionType>::dir2>(modStr,rem,possibleToBePossitive,vals);
		reasoninInDirection<RemovalSquare,ModuleSquare,-Direction<DirectionType>::dir3>(modStr,rem,possibleToBePossitive,vals);
		reasoninInDirection<RemovalSquare,ModuleSquare,-Direction<DirectionType>::dir4>(modStr,rem,possibleToBePossitive,vals);
	}
};

///end HysteresisThresholding

template <typename ValueType,typename CoordinateType,int Size>
	ValueType gaussian5x5(CoordinateType x,CoordinateType y)
	{
		static const ValueType tab[Size][Size]={{2,4,5,4,2},{2,9,12,9,4},{5,12,15,12,2},{4,9,12,9,4},{2,4,5,4,2}};
		return tab[y+2][x+2];
	}
template <typename ValueType,typename CoordinateType,int Size>
	ValueType gaussianLoG(CoordinateType x,CoordinateType y)
	{
		static const ValueType tab[Size][Size]={{0,0,-1,0,0},{0,-1,-2,-1,0},{-1,-2,16,-2,-1},{0,-1,-2,-1,0},{0,0,-1,0,0}};
		//static const ValueType tab[Size][Size]={{0,0,1,0,0},{0,1,2,1,0},{1,2,-16,2,1},{0,1,2,1,0},{0,0,1,0,0}};
		return tab[y+2][x+2];
	}

class CannyEdgeDetection
{
public:
	

	template <typename IM,class ValueType>
	static void detectEdges(IM& im,ValueType low,ValueType high)
	{
		//Typy pomocnicze
		static const int Size=2;
		typedef ValueType GradientStrengthType;
		typedef char GradientDirectionType;
		typedef int SumAccumlatorType;
		typedef IM::ValueType PixelType;
		typedef AllocatedValueMap<GradientDirectionType,IM::CoordinateType> GradMap;
		typedef AllocatedValueMap<GradientStrengthType,IM::CoordinateType> StMap;
		typedef typename IM::Iterator<Size>::SquareType SquareType;
		typedef AllocatedValueMap<bool,IM::CoordinateType> BoolMap;
		typename typedef numeric_limits<PixelType> limit;
		static const int NSS=3;//NEIGHBORHOOD_SQUARE_SIZE=3
		static const int GMS=5;//GAUSSIAN_MASK_SIZE=5 
		

		HLvalues<GradientStrengthType> values;
		values.high=high;
		values.low=low;

	  typename typedef Mask<GMS,IM::ValueType,IM::CoordinateType,&gaussian5x5<IM::ValueType,IM::CoordinateType,GMS > > GMask;

	  {//filtr Gaussa iterujemy po oknach 5x5
			for(IM::Iterator<GMS> iter=im.begin<GMS>();iter!=im.end<GMS>();iter++)
				filter<GMask,IM::Iterator<GMS>::SquareType,SumAccumlatorType>(*iter);
	  }//filtr Gaussa

	  //Mapy wartoœci
	  
	  StMap gradientStrengthMap(im);//modu³y, si³a gradientu
	  GradMap edgeDirectionMap(im);//kierunek krawêdzi, do której ewnetualnie nale¿y piksel
	  BoolMap pixelsToRemove(im);//mapa oznaczeñ pikseli do usuniêcia - wartoœæ true to usun¹æ

	  edgeDirectionMap.init();//
	  pixelsToRemove.init();
{
GradMap gradientDirectionMap(im);//kierunek gradientów
	  {////////////////////////// okna 2x2
						GradMap::Iterator<Size> gradientDirectionIter=gradientDirectionMap.begin<Size>();
						StMap::Iterator<Size> gradientStrengthIter=gradientStrengthMap.begin<Size>();
						for(IM::Iterator<Size> iter=im.begin<Size>();iter!=im.end<Size>();iter++)
						{
							Derivative<GradientStrengthType> der=DerivativePolicy<SquareType,GradientStrengthType>::getDerivative((*iter));
							(*gradientDirectionIter)(0,0)=
							NonMaximumSuppressionAlgorithms::direction<GradientStrengthType,GradientDirectionType>(atan(der.Gy/der.Gx));
							(*gradientStrengthIter)(0,0)=sqrt((der.Gx*der.Gx)+(der.Gy*der.Gy));
							gradientDirectionIter++;
							gradientStrengthIter++;
						}
	  }//////////////////////////

		{/////////////////////////// okna 3x3
						GradMap::Iterator<NSS> sIter=gradientDirectionMap.begin<NSS>();
						GradMap::Iterator<NSS> dIter=edgeDirectionMap.begin<NSS>();
						/*Sprawdzamy czy przez piksel przechodzi krawêdŸ*/
								for(;sIter!=gradientDirectionMap.end<NSS>();sIter++)
								{
										(*dIter)(0,0)=NonMaximumSuppressionAlgorithms::centerPixelBuildsEdge(*sIter);
										dIter++;
								}
		}/////////////////////////////
}


		{///////////////////t³umienie niemaksymalne okna 3x3
					
					StMap::Iterator<NSS> gradientStrengthIter=gradientStrengthMap.begin<NSS>();//
					BoolMap::Iterator<NSS> removalIter=pixelsToRemove.begin<NSS>();//
					GradMap::Iterator<NSS> edgeIter=edgeDirectionMap.begin<NSS>();
					for(;edgeIter!=edgeDirectionMap.end<NSS>();edgeIter++)
					{
						if((*edgeIter)(0,0)!=Direction<GradientDirectionType>::none)
						{
							(*removalIter)(0,0)=NonMaximumSuppressionAlgorithms::classifiesToRemoval(*gradientStrengthIter,(*edgeIter)(0,0));
						}else
						{
							(*removalIter)(0,0)=true;
						}
						gradientStrengthIter++;
						removalIter++;
					}
					
		}////////////////////////
	
		//progowanie dla high - okna 3x3
		BoolMap possibleCandidates(pixelsToRemove,true);//kopia 
		{
					StMap::Iterator<NSS> gradientStrengthIter=gradientStrengthMap.begin<NSS>();//ok
					BoolMap::Iterator<NSS> removalIter=pixelsToRemove.begin<NSS>();//ok
					for(;gradientStrengthIter!=gradientStrengthMap.end<NSS>();gradientStrengthIter++)
					{
						if((*gradientStrengthIter)(0,0)<values.high)
							(*removalIter)(0,0)=true;
						removalIter++;
					}
		}
		

		{// Hysteresis Thresholding
			//okna 3x3
				BoolMap firstReasonedEdges(pixelsToRemove,true);
				BoolMap::Iterator<NSS> destRemovalsIter=pixelsToRemove.begin<NSS>();
				BoolMap::Iterator<NSS> possibleEdgesIter=possibleCandidates.begin<NSS>();
				StMap::Iterator<NSS> gradientStrengthIter=gradientStrengthMap.begin<NSS>();//ok
				for(BoolMap::Iterator<NSS> firstRemoval=firstReasonedEdges.begin<NSS>();firstRemoval!=firstReasonedEdges.end<NSS>();firstRemoval++)
				{
					if(!(*firstRemoval)(0,0))
					{
						HysteresisThresholding<GradientDirectionType>::reasonPosotiveNeighbors
							(*gradientStrengthIter,*destRemovalsIter,*possibleEdgesIter,values);   
					}
					gradientStrengthIter++;
					destRemovalsIter++;
					possibleEdgesIter++;
				}
				
		}

		{//fina³
				BoolMap::Iterator<1> bmRemoval=pixelsToRemove.begin<1>();
				for(IM::Iterator<1> iter=im.begin<1>();iter!=im.end<1>();iter++)
				{

					if((*bmRemoval)(0,0)) (*iter)(0,0)=PixelType();
					else
						#undef max
						(*iter)(0,0)=limit::max();
						#define max
					bmRemoval++;
				}
		}
	}
};

class LaplacianOfGaussian
{
public:
	template <typename IM,class ValueType>
	static void apply(IM& im,ValueType low,ValueType high)
	{
		static const int Size=2;
		typedef ValueType GradientStrengthType;
		typedef char GradientDirectionType;
		typedef int SumAccumlatorType;
		typedef IM::ValueType PixelType;
		
		typename typedef numeric_limits<PixelType> limit;
		static const int NSS=3;//NEIGHBORHOOD_SQUARE_SIZE=3
		static const int ZEROCROSSING=2;
		static const int GMS=5;//GAUSSIAN_MASK_SIZE=5 
		typedef AllocatedValueMap<int,IM::CoordinateType> ZeroCrossingMapType;

		typedef AllocatedValueMap<bool,IM::CoordinateType> BoolMap;
	    typename typedef Mask<GMS,IM::ValueType,IM::CoordinateType,&gaussian5x5<IM::ValueType,IM::CoordinateType,GMS > > GMask;
		typename typedef Mask<GMS,SumAccumlatorType,IM::CoordinateType,&gaussianLoG<SumAccumlatorType,IM::CoordinateType,GMS > > LoGMask;

	  
	  {//filtr Gaussa iterujemy po oknach 5x5
			for(IM::Iterator<GMS> iter=im.begin<GMS>();iter!=im.end<GMS>();iter++)
				filter<GMask,IM::Iterator<GMS>::SquareType,SumAccumlatorType>(*iter);
	  }//filtr Gaussa
		ZeroCrossingMapType zeroCrossingMap(im,true);

		ZeroCrossingMapType zeroCrossingMapSource(im,true);

		BoolMap crossedZeroValueMap(im);
		crossedZeroValueMap.init();
			  {
				ZeroCrossingMapType::Iterator<GMS> iterSource=zeroCrossingMapSource.begin<GMS>();
					for(ZeroCrossingMapType::Iterator<GMS> iter=zeroCrossingMap.begin<GMS>();iter!=zeroCrossingMap.end<GMS>();iter++)
					{
						filter_log<LoGMask,ZeroCrossingMapType::Iterator<GMS>::SquareType,SumAccumlatorType>(*iterSource,*iter);
						iterSource++;
					}
			  }
			{
				BoolMap::Iterator<ZEROCROSSING> boolIterator=crossedZeroValueMap.begin<ZEROCROSSING>();
				for(ZeroCrossingMapType::Iterator<ZEROCROSSING> iter=zeroCrossingMap.begin<ZEROCROSSING>();iter!=zeroCrossingMap.end<ZEROCROSSING>();iter++)
				{
					(*boolIterator)(0,0)=ZeroCrossing::crossesZero2x2(*iter,low,high);
					boolIterator++;
				}
			}
			{
				/*BoolMap::Iterator<2> boolIterator=crossedZeroValueMap.begin<2>();
				for(ZeroCrossingMapType::Iterator<2> iter=zeroCrossingMap.begin<2>();iter!=zeroCrossingMap.end<2>();iter++)
				{
					(*boolIterator)(0,0)=ZeroCrossing::crossesZero2x2(*iter,low,high)&&(*boolIterator)(0,0);
					boolIterator++;
				}*/
			}
		{//fina³
				BoolMap::Iterator<1> bmRemoval=crossedZeroValueMap.begin<1>();
				for(IM::Iterator<1> iter=im.begin<1>();iter!=im.end<1>();iter++)
				{

					if(!(*bmRemoval)(0,0)) (*iter)(0,0)=PixelType();
					else
						#undef max
						(*iter)(0,0)=limit::max();
						#define max
					bmRemoval++;
				}
		}
	}

};
