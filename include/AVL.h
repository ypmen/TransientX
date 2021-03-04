/*
 * AVL.h
 *
 *  Created on: 2018年3月14日
 *      Author: ypmen
 */

#ifndef AVL_H_
#define AVL_H_

#include <iostream>
#include <math.h>
using namespace std;

template <typename T> class AVLNode
{
public:
	AVLNode();
	AVLNode(T val);
	AVLNode(T value, AVLNode<T> * left, AVLNode<T> * right);
	~AVLNode();
	void setleft(AVLNode<T> * node);
	void setright(AVLNode<T> * node);
	AVLNode<T> * getleft(void){return left;}
	AVLNode<T> * getright(void){return right;}
	int getheight(void){return height;}
	int getcount(void){return count;}
public:
	T value;
private:
	int height;
	int count;
	AVLNode<T> * left;
	AVLNode<T> * right;
};

template <typename T> class AVLTree
{
public:
	AVLTree();
	AVLTree(AVLNode<T> * node);
	AVLTree(T * seq, long int size);
	~AVLTree();
	void insertValue(T val);
	void removeValue(T val);
	T minValue(void);
	T maxValue(void);
	void printValues(void);
	void setCountEqual(void);
	T getMedian(void);
private:
	AVLNode<T> * newNode(T val);
	AVLNode<T> * insertNode(AVLNode<T> * node, T val);
	AVLNode<T> * removeNode(AVLNode<T> * node, T val);
	AVLNode<T> * minValueNode(AVLNode<T> * node);
	AVLNode<T> * maxValueNode(AVLNode<T> * node);
	AVLNode<T> * rightRotate(AVLNode<T> * node);
	AVLNode<T> * leftRotate(AVLNode<T> * node);
	int getBalance(AVLNode<T> * node);
	int getCountCompare(AVLNode<T> * node);
	void preOrder(AVLNode<T> * node);
private:
	AVLNode<T> * root;
};

#endif /* AVL_H_ */
