#pragma once

#include "AADTape.h"
#include <cmath>

//  Base CRTP class so operators catch expressions
template <class E>
struct Expression
{
    double value() const { return static_cast<const E*>(this)->value(); }
};

//  Note that Number is a leaf expression
//  Defined in the bottom of the file

//  Binary expression
//  LHS : the expression on the left
//  RHS : the expression on the right
//  OP : the binary function 
template <class LHS, class RHS, class OP>
class BinaryExpression 
    //  CRTP
    : public Expression<BinaryExpression<LHS, RHS, OP>>
{
    const double myValue;

    const LHS lhs;
    const RHS rhs;

public:

    //  Constructor out of 2 expressions
    //  Note : eager evaluation on construction
    explicit BinaryExpression(
        const Expression<LHS>& l,
        const Expression<RHS>& r)
        : myValue(OP::eval(l.value(), r.value())), 
        lhs(static_cast<const LHS&>(l)), 
        rhs(static_cast<const RHS&>(r))
    {}

    //  Value accessors
    double value() const { return myValue; }

    //	Expression template magic
    //  Expressions know
    //  AT COMPILE TIME
    //  the number of numbers in their sub-DAG
    enum { numNumbers = LHS::numNumbers + RHS::numNumbers };

    //  Push adjoint down the DAG
    //  And write information into the node for the expression
    //  N : total number of numbers in the expression
    //  n : numbers already processed
    template <size_t N, size_t n>
    void pushAdjoint(
        //  Node for the complete expression being processed
        MultiNode<N>& exprNode,     
        //  Adjoint cumulated for this binary node
        const double adjoint)       
        const
    {
        //  Push on the left, if numbers there
        if (LHS::numNumbers > 0)
        {
            lhs.pushAdjoint<N, n>(
                exprNode, 
                adjoint * OP::leftDerivative(lhs.value(), rhs.value(), value()));
        }

        //  Push on the right
        if (RHS::numNumbers > 0)
        {
            //  Note left push processed LHS::numNumbers numbers
            //  So the next number to be processed is n + LHS::numNumbers
            rhs.pushAdjoint<N, n + LHS::numNumbers>(
                exprNode, 
                adjoint * OP::rightDerivative(lhs.value(), rhs.value(), value()));
        }
    }
};

//  "Concrete" binaries, we only need to define operations and derivatives
struct OPMult
{
    static const double eval(const double l, const double r) 
    { 
        return l * r; 
    }
    
    static const double leftDerivative
        (const double l, const double r, const double v) 
    { 
        return r; 
    }
    
    static const double rightDerivative
        (const double l, const double r, const double v) 
    { 
        return l; 
    }
};

struct OPAdd
{
    static const double eval(const double l, const double r)
    { 
        return l + r; 
    }
    
    static const double leftDerivative
        (const double l, const double r, const double v)
    { 
        return 1.0; 
    }
    
    static const double rightDerivative
        (const double l, const double r, const double v)
    { 
        return 1.0; 
    }
};

struct OPSub
{
    static const double eval(const double l, const double r)
    {
        return l - r;
    }

    static const double leftDerivative
    (const double l, const double r, const double v)
    {
        return 1.0;
    }

    static const double rightDerivative
    (const double l, const double r, const double v)
    {
        return -1.0;
    }
};

struct OPDiv
{
    static const double eval(const double l, const double r)
    {
        return l / r;
    }

    static const double leftDerivative
    (const double l, const double r, const double v)
    {
        return 1.0 / r;
    }

    static const double rightDerivative
    (const double l, const double r, const double v)
    {
        return -l / r / r;
    }
};

struct OPPow
{
    static const double eval(const double l, const double r)
    {
        return pow(l, r);
    }

    static const double leftDerivative
    (const double l, const double r, const double v)
    {
        return r*v / l;
    }

    static const double rightDerivative
    (const double l, const double r, const double v)
    {
        return log(l)*v;
    }
};

//  Operator overloading for binary expressions
//  So DAG is built on the stack at compile time
//  And traversed at run time for evaluation and propagation

template <class LHS, class RHS>
 BinaryExpression<LHS, RHS, OPMult> operator*(
    const Expression<LHS>& lhs, const Expression<RHS>& rhs) 
{ 
    return BinaryExpression<LHS, RHS, OPMult>(lhs, rhs); 
}

template <class LHS, class RHS>
 BinaryExpression<LHS, RHS, OPAdd> operator+(
    const Expression<LHS>& lhs, const Expression<RHS>& rhs)
{
    return BinaryExpression<LHS, RHS, OPAdd>(lhs, rhs);
}

template <class LHS, class RHS>
 BinaryExpression<LHS, RHS, OPSub> operator-(
    const Expression<LHS>& lhs, const Expression<RHS>& rhs)
{
    return BinaryExpression<LHS, RHS, OPSub>(lhs, rhs);
}

template <class LHS, class RHS>
 BinaryExpression<LHS, RHS, OPDiv> operator/(
    const Expression<LHS>& lhs, const Expression<RHS>& rhs)
{
    return BinaryExpression<LHS, RHS, OPDiv>(lhs, rhs);
}

template <class LHS, class RHS>
 BinaryExpression<LHS, RHS, OPPow> pow(
    const Expression<LHS>& lhs, const Expression<RHS>& rhs)
{
    return BinaryExpression<LHS, RHS, OPPow>(lhs, rhs);
}

//  Unary expressions : Same logic with one argument

//  The CRTP class
template <class ARG, class OP>
class UnaryExpression
    //  CRTP
    : public Expression<UnaryExpression<ARG, OP>>
{
    const double myValue;

    const ARG arg;
    //  For binary operators with a double on one side
    const double dArg = 0.0;   

public:

    //  Constructor 
    //  Note : eager evaluation on construction
    explicit UnaryExpression(
        const Expression<ARG>& a)
        : myValue(OP::eval(a.value(), 0.0)), arg(static_cast<const ARG&>(a)) {}

    //  Constructor for binary expressions with a double on one side
    explicit UnaryExpression(
        const Expression<ARG>& a,
        const double b)
        : myValue(OP::eval(a.value(), b)), arg(static_cast<const ARG&>(a)), dArg(b) {}

    //  Value accessors
    double value() const { return myValue; }

    //	Expression template magic
    enum { numNumbers = ARG::numNumbers };

    //  Push adjoint down the expression DAG
    template <size_t N, size_t n>
    void pushAdjoint(
        MultiNode<N>& exprNode,     //  Node for the complete expression
        const double adjoint)       //  Adjoint cumulated on the node
        const
    {
        //  Push to argument, if numbers there
        if (ARG::numNumbers > 0)
        {
            arg.pushAdjoint<N, n>(
                exprNode, 
                adjoint * OP::derivative(arg.value(), value(), dArg));
        }
    }
};

//  The unary operators

struct OPExp
{
    static const double eval(const double r, const double d) 
    { 
        return exp(r); 
    }
    
    static const double derivative
        (const double r, const double v, const double d)
    { 
        return v; 
    }
};

struct OPLog
{
    static const double eval(const double r, const double d)
    { 
        return log(r); 
    }
    
    static const double derivative
        (const double r, const double v, const double d)
    { 
        return 1.0 / r; 
    }
};

struct OPSqrt
{
    static const double eval(const double r, const double d)
    { 
        return sqrt(r); 
    }

    static const double derivative
        (const double r, const double v, const double d)
    { 
        return 0.5 / v; 
    }
};

//  Binary operators with a double on one side 

//  * double or double *
struct OPMultD
{
    static const double eval(const double r, const double d)
    {
        return r * d;
    }

    static const double derivative
    (const double r, const double v, const double d)
    {
        return d;
    }
};

//  + double or double +
struct OPAddD
{
    static const double eval(const double r, const double d)
    {
        return r + d;
    }

    static const double derivative
    (const double r, const double v, const double d)
    {
        return 1.0;
    }
};

//  double -
struct OPSubDL
{
    static const double eval(const double r, const double d)
    {
        return d - r;
    }

    static const double derivative
    (const double r, const double v, const double d)
    {
        return -1.0;
    }
};

//  - double
struct OPSubDR
{
    static const double eval(const double r, const double d)
    {
        return r - d;
    }

    static const double derivative
    (const double r, const double v, const double d)
    {
        return 1.0;
    }
};

//  double /
struct OPDivDL
{
    static const double eval(const double r, const double d)
    {
        return d / r;
    }

    static const double derivative
    (const double r, const double v, const double d)
    {
        return -d / r / r;
    }
};

//  / double
struct OPDivDR
{
    static const double eval(const double r, const double d)
    {
        return r / d;
    }

    static const double derivative
    (const double r, const double v, const double d)
    {
        return 1.0 / d;
    }
};

//  pow (d/)
struct OPPowDL
{
    static const double eval(const double r, const double d)
    {
        return pow(d, r);
    }

    static const double derivative
    (const double r, const double v, const double d)
    {
        return log(d) * v;
    }
};

//  pow (/d)
struct OPPowDR
{
    static const double eval(const double r, const double d)
    {
        return pow(r, d);
    }

    static const double derivative
    (const double r, const double v, const double d)
    {
        return d * v / r;
    }
};

//  And overloading

template <class ARG>
 UnaryExpression<ARG, OPExp> exp(const Expression<ARG>& arg)
{ 
    return UnaryExpression<ARG, OPExp>(arg);
}

template <class ARG>
 UnaryExpression<ARG, OPLog> log(const Expression<ARG>& arg)
{
    return UnaryExpression<ARG, OPLog>(arg);
}

template <class ARG>
 UnaryExpression<ARG, OPSqrt> sqrt(const Expression<ARG>& arg)
{
    return UnaryExpression<ARG, OPSqrt>(arg);
}

//  Binary operators with a double on one side 

template <class ARG>
 UnaryExpression<ARG, OPMultD> operator*(
    const double d, const Expression<ARG>& rhs)
{
    return UnaryExpression<ARG, OPMultD>(rhs, d);
}

template <class ARG>
 UnaryExpression<ARG, OPMultD> operator*(
    const Expression<ARG>& lhs, const double d)
{
    return UnaryExpression<ARG, OPMultD>(lhs, d);
}

template <class ARG>
 UnaryExpression<ARG, OPAddD> operator+(
    const double d, const Expression<ARG>& rhs)
{
    return UnaryExpression<ARG, OPAddD>(rhs, d);
}

template <class ARG>
 UnaryExpression<ARG, OPAddD> operator+(
    const Expression<ARG>& lhs, const double d)
{
    return UnaryExpression<ARG, OPAddD>(lhs, d);
}

template <class ARG>
 UnaryExpression<ARG, OPSubDL> operator-(
    const double d, const Expression<ARG>& rhs)
{
    return UnaryExpression<ARG, OPSubDL>(rhs, d);
}

template <class ARG>
 UnaryExpression<ARG, OPSubDR> operator-(
    const Expression<ARG>& lhs, const double d)
{
    return UnaryExpression<ARG, OPSubDR>(lhs, d);
}

template <class ARG>
 UnaryExpression<ARG, OPDivDL> operator/(
    const double d, const Expression<ARG>& rhs)
{
    return UnaryExpression<ARG, OPDivDL>(rhs, d);
}

template <class ARG>
 UnaryExpression<ARG, OPDivDR> operator/(
    const Expression<ARG>& lhs, const double d)
{
    return UnaryExpression<ARG, OPDivDR>(lhs, d);
}

template <class ARG>
 UnaryExpression<ARG, OPPowDL> pow(
    const double d, const Expression<ARG>& rhs)
{
    return UnaryExpression<ARG, OPPowDL>(rhs, d);
}

template <class ARG>
 UnaryExpression<ARG, OPPowDR> pow(
    const Expression<ARG>& lhs, const double d)
{
    return UnaryExpression<ARG, OPPowDR>(lhs, d);
}

//  Comparison, as normal

template<class E, class F>
 bool operator==(const Expression<E>& lhs, const Expression<F>& rhs)
{
    return lhs.value() == rhs.value();
}
template<class E>
 bool operator==(const Expression<E>& lhs, const double& rhs)
{
    return lhs.value() == rhs;
}
template<class E>
 bool operator==(const double& lhs, const Expression<E>& rhs)
{
    return lhs == rhs.value();
}

template<class E, class F>
 bool operator!=(const Expression<E>& lhs, const Expression<F>& rhs)
{
    return lhs.value() != rhs.value();
}
template<class E>
 bool operator!=(const Expression<E>& lhs, const double& rhs)
{
    return lhs.value() != rhs;
}
template<class E>
 bool operator!=(const double& lhs, const Expression<E>& rhs)
{
    return lhs != rhs.value();
}

template<class E, class F>
 bool operator<(const Expression<E>& lhs, const Expression<F>& rhs)
{
    return lhs.value() < rhs.value();
}
template<class E>
 bool operator<(const Expression<E>& lhs, const double& rhs)
{
    return lhs.value() < rhs;
}
template<class E>
 bool operator<(const double& lhs, const Expression<E>& rhs)
{
    return lhs < rhs.value();
}

template<class E, class F>
 bool operator>(const Expression<E>& lhs, const Expression<F>& rhs)
{
    return lhs.value() > rhs.value();
}
template<class E>
 bool operator>(const Expression<E>& lhs, const double& rhs)
{
    return lhs.value() > rhs;
}
template<class E>
 bool operator>(const double& lhs, const Expression<E>& rhs)
{
    return lhs > rhs.value();
}

template<class E, class F>
 bool operator<=(const Expression<E>& lhs, const Expression<F>& rhs)
{
    return lhs.value() <= rhs.value();
}
template<class E>
 bool operator<=(const Expression<E>& lhs, const double& rhs)
{
    return lhs.value() <= rhs;
}
template<class E>
 bool operator<=(const double& lhs, const Expression<E>& rhs)
{
    return lhs <= rhs.value();
}

template<class E, class F>
 bool operator>=(const Expression<E>& lhs, const Expression<F>& rhs)
{
    return lhs.value() >= rhs.value();
}
template<class E>
 bool operator>=(const Expression<E>& lhs, const double& rhs)
{
    return lhs.value() >= rhs;
}
template<class E>
 bool operator>=(const double& lhs, const Expression<E>& rhs)
{
    return lhs >= rhs.value();
}

//  Finally, unary +/- operators

template <class RHS>
UnaryExpression<RHS, OPSubDL> operator-
(const Expression<RHS>& rhs)
{
    return 0.0 - rhs;
}

template <class RHS>
Expression<RHS> operator+
(const Expression<RHS>& rhs)
{
    return rhs;
}

//  The Number type, also an expression

class Number : public Expression<Number>
{
    //  The value and node for this number, as normal
    double myValue;
    Node* myNode;

    //  Node creation on tape, as normal
    template <size_t N>
    MultiNode<N>* createMultiNode()
    {
        //  Placement syntax to allocate in place on tape
        return new (tape->allocate(sizeof(MultiNode<N>))) MultiNode<N>;
    }

    //  This is where, on assignment or construction from an expression,
    //      that derivatives are pushed through the expression's DAG 
    //      into the node
    template<class E>
    void fromExpr(const Expression<E>& e)
    {
        //  Build node
        auto* node = createMultiNode<E::numNumbers>();
        
        //  Push adjoints through expression DAG from 1 on top
        static_cast<const E&>(e).pushAdjoint<E::numNumbers, 0>(*node, 1.0);

        //  Set my node
        myNode = node;
    }

public:

    //  Expression template magic
    enum { numNumbers = 1 };

    //  Push adjoint
    //  This is where the derivatives and pointers are set on the node
    //  Push adjoint down the expression DAG
    template <size_t N, size_t n>
    void pushAdjoint(
        MultiNode<N>& exprNode,     //  Node for the complete expression
        const double adjoint)       //  Adjoint cumulated on the node
        const
    {
        exprNode.arguments[n] = myNode;
        exprNode.derivatives[n] = adjoint;
    }

    //  Static access to tape, as normal
    static thread_local Tape* tape;

    //  Constructors

    Number() {}

    explicit Number(const double val) : myValue(val)
    {
        //  Create leaf
        myNode = createMultiNode<0>();
    }

    Number& operator=(const double val)
    {
        myValue = val;
        myNode = createMultiNode<0>();
        return *this;
    }

    //  No need for copy and assignment
    //  Default ones do the right thing:
    //      copy value and pointer to node

    //  Construct or assign from expression
    
    template <class E>
    Number(const Expression<E>& e) : myValue(e.value())
    {
        fromExpr<E>(static_cast<const E&>(e));
    }

    template <class E>
    Number& operator=
        (const Expression<E>& e) 
    {
        myValue = e.value();
        fromExpr<E>(static_cast<const E&>(e));
        return *this;
    }

    //  All the normal accessors and propagators, as in normal AAD code
    
    //  Put on tape
    void putOnTape()
    {
        myNode = createMultiNode<0>();
    }

    //  Accessors: value and adjoint

    double& value()
    {
        return myValue;
    }
    double value() const
    {
        return myValue;
    }
    double& adjoint()
    {
        return myNode->adjoint;
    }
    double adjoint() const
    {
        return myNode->adjoint;
    }

    //  Propagation

    //  Reset all adjoints on the tape
    static void resetAdjoints()
    {
        for (Node& node : *tape) node.adjoint = 0.0;
    }
    //  Propagate adjoints
    //      from and to both INCLUSIVE
    static void propagateAdjoints(
        Tape::iterator propagateFrom,
        Tape::iterator propagateTo)
    {
        auto it = propagateFrom;
        while (it != propagateTo)
        {
            it->propagate();
            --it;
        }
        it->propagate();
    }

    //  Convenient overloads

    //  Set the adjoint on this node to 1,
    //  Then propagate from the node
    void propagateAdjoints(
        //  We start on this number's node
        Tape::iterator propagateTo,
        //  reset adjoints first?
        const bool reset = false)
    {
        //  Reset
        if (reset) resetAdjoints();
        //  Set this adjoint to 1
        adjoint() = 1.0;
        //  Find node on tape
        auto it = tape->find(myNode);
        //  Reverse and propagate until we hit the stop
        while (it != propagateTo)
        {
            it->propagate();
            --it;
        }
        it->propagate();
    }

    //  These 2 set the adjoint to 1 on this node
    void propagateToStart(
        const bool reset = false)
    {
        propagateAdjoints(tape->begin(), reset);
    }
    void propagateToMark(
        const bool reset = false)
    {
        propagateAdjoints(tape->markIt(), reset);
    }

    //  This one only propagates
    //  Note: propagation starts at mark - 1
    static void propagateMarkToStart()
    {
        propagateAdjoints(--tape->markIt(), tape->begin());
    }

    //  Unary operators

    template <class E>
    Number& operator+=(const Expression<E>& e)
    {
        *this = *this + e;
        return *this;
    }

    template <class E>
    Number& operator*=(const Expression<E>& e)
    {
        *this = *this * e;
        return *this;
    }

    template <class E>
    Number& operator-=(const Expression<E>& e)
    {
        *this = *this - e;
        return *this;
    }

    template <class E>
    Number& operator/=(const Expression<E>& e)
    {
        *this = *this / e;
        return *this;
    }
};

