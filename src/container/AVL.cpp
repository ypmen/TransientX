//============================================================================
// Name        : AVL.cpp
// Author      : ypmen
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include "AVL.h"
#include <cmath>

template <typename T>
AVLTree<T>::AVLTree()
{
	root = NULL;
}

template <typename T>
AVLTree<T>::AVLTree(AVLNode<T> * node)
{
	root = node;
}

template <typename T>
AVLTree<T>::AVLTree(T * seq, long int size)
{
	root = NULL;
	for (long int i=0; i< size; i++)
	{
		insertValue(seq[i]);
	}
}

template <typename T>
AVLTree<T>::~AVLTree()
{
	if (root != NULL)
	{
		 delete root;
		 root = NULL;
	}
}

template <typename T>
void AVLTree<T>::insertValue(T val)
{
	root = insertNode(root, val);
	//cout<<endl;
}

template <typename T>
AVLNode<T> * AVLTree<T>::insertNode(AVLNode<T> * node, T val)
{
	if (node == NULL)
		return newNode(val);

	if (val <= node->value)
	{
		//cout<<"l ";
		node->setleft(insertNode(node->getleft(), val));
	}
	else if (val > node->value)
	{
		//cout<<"r ";
		node->setright(insertNode(node->getright(), val));
	}
	/*else
	{
		//cout<<"= ";
		AVLNode<T> * tmp = newNode(val);
		tmp->setleft(node->getleft());
		node->setleft(tmp);
	}*/

	int balance = getBalance(node);

	AVLNode<T> * tmpleft = node->getleft();
	AVLNode<T> * tmpright = node->getright();

	//Left Left Case
	if (balance > 1 and val <= tmpleft->value)
		return rightRotate(node);
	//Right Right Case
	if (balance < -1 and val > tmpright->value)
		return leftRotate(node);
	//Left Right Case
	if (balance > 1 and val > tmpleft->value)
	{
		node->setleft(leftRotate(tmpleft));
		return rightRotate(node);
	}
	//Right Left Case
	if (balance < -1 and val <= tmpright->value)
	{
		node->setright(rightRotate(tmpright));
		return leftRotate(node);
	}

	return node;
}

template <typename T>
AVLNode<T> * AVLTree<T>::newNode(T val)
{
	AVLNode<T> * node = new AVLNode<T>(val);
	return node;
}


template <typename T>
void AVLTree<T>::removeValue(T val)
{
	root = removeNode(root, val);
	//cout<<endl;
}

template <typename T>
AVLNode<T> * AVLTree<T>::removeNode(AVLNode<T> * node, T val)
{
	if (node == NULL)
		return node;

	if (val < node->value)
	{
		//cout<<"l ";
		node->setleft(removeNode(node->getleft(), val));
	}
	else if (val > node->value)
	{
		//cout<<"r ";
		node->setright(removeNode(node->getright(), val));
	}
	else
	{
		//cout<<"= ";
		//cout<<endl;
		if (node->getleft() == NULL or node->getright() == NULL)
		{
			AVLNode<T> * tmp = node->getleft() ? node->getleft() : node->getright();
			if (tmp == NULL)
			{
				tmp = node;
				node = NULL;
			}
			else
			{
				AVLNode<T> * tmp2 = node;
				node = tmp;
				tmp = tmp2;
			}
			tmp->setleft(NULL);
			tmp->setright(NULL);
			delete tmp;
		}
		else
		{
			if (getCountCompare(node) <= 0)
			{
				AVLNode<T> * tmp = minValueNode(node->getright());
				node->value = tmp->value;
				node->setright(removeNode(node->getright(), tmp->value));
			}
			else
			{
				AVLNode<T> * tmp = maxValueNode(node->getleft());
				node->value = tmp->value;
				node->setleft(removeNode(node->getleft(), tmp->value));
			}
		}
	}

	if (node == NULL)
		return node;

	int balance = getBalance(node);

	AVLNode<T> * tmpleft = node->getleft();
	AVLNode<T> * tmpright = node->getright();

	//Left Left Case
	if (balance > 1 and getBalance(node->getleft()) >= 0)
		return rightRotate(node);
	//Right Right Case
	if (balance < -1 and getBalance(node->getright()) <= 0)
		return leftRotate(node);
	//Left Right Case
	if (balance > 1 and getBalance(node->getleft()) < 0)
	{
		node->setleft(leftRotate(tmpleft));
		return rightRotate(node);
	}
	//Right Left Case
	if (balance < -1 and getBalance(node->getright()) > 0)
	{
		node->setright(rightRotate(tmpright));
		return leftRotate(node);
	}

	return node;
}

template <typename T>
AVLNode<T> * AVLTree<T>::minValueNode(AVLNode<T> * node)
{
	AVLNode<T> * curr = node;

	while (curr->getleft() != NULL)
		curr = curr->getleft();

	return curr;
}

template <typename T>
AVLNode<T> * AVLTree<T>::maxValueNode(AVLNode<T> * node)
{
	AVLNode<T> * curr = node;

	while (curr->getright() != NULL)
		curr = curr->getright();

	return curr;
}

template <typename T>
T AVLTree<T>::minValue(void)
{
	AVLNode<T> * minNode = minValueNode(root);
	return minNode->value;
}

template <typename T>
T AVLTree<T>::maxValue(void)
{
	AVLNode<T> * maxNode = maxValueNode(root);
	return maxNode->value;
}

template <typename T>
void AVLTree<T>::printValues(void)
{
	preOrder(root);
	cout<<endl;
}

template <typename T>
void AVLTree<T>::preOrder(AVLNode<T> * node)
{
	if(node != NULL)
	{
		cout<<node->value<<" ";
		//cout<<node->getheight()<<" ";
		//cout<<node->getcount()<<" ";
		//cout<<endl;
		preOrder(node->getleft());
		preOrder(node->getright());
	}
}

template <typename T>
AVLNode<T> * AVLTree<T>::rightRotate(AVLNode<T> * node)
{
	AVLNode<T> * x = node->getleft();
	AVLNode<T> * y = x->getright();

	node->setleft(y);
	x->setright(node);

	return x;
}

template <typename T>
AVLNode<T> * AVLTree<T>::leftRotate(AVLNode<T> * node)
{
	AVLNode<T> * x = node->getright();
	AVLNode<T> * y = x->getleft();

	node->setright(y);
	x->setleft(node);

	return x;
}

template <typename T>
int AVLTree<T>::getBalance(AVLNode<T> * node)
{
	if (node == NULL)
		return 0;
	AVLNode<T> * tmpleft = node->getleft();
	AVLNode<T> * tmpright = node->getright();
	int lhei = tmpleft ? tmpleft->getheight() : -1;
	int rhei = tmpright ? tmpright->getheight() : -1;
	return lhei-rhei;
}

template <typename T>
int AVLTree<T>::getCountCompare(AVLNode<T> * node)
{
	if (node == NULL)
		return 0;

	AVLNode<T> * tmpleft = node->getleft();
	AVLNode<T> * tmpright = node->getright();

	int lcon = 0;
	int rcon = 0;
	if (tmpleft != NULL)
		lcon = tmpleft->getcount();
	if (tmpright != NULL)
		rcon = tmpright->getcount();

	return lcon - rcon;
}

template <typename T>
void AVLTree<T>::setCountEqual(void)
{
	int com = 0;
	while (abs(com = getCountCompare(root)) > 1)
	{
		T tmp = root->value;
		root = removeNode(root, tmp);
		if (com > 1)
			root->setright(insertNode(root->getright(), tmp));
		else
			root->setleft(insertNode(root->getleft(), tmp));
	}
}

template <typename T>
T AVLTree<T>::getMedian(void)
{
	setCountEqual();
	int com = getCountCompare(root);
	AVLNode<T> * tmp = NULL;

	if (com == 1)
	{
		tmp = maxValueNode(root->getleft());
		return (root->value + tmp->value)/2;
	}
	else if (com == -1)
	{
		tmp = minValueNode(root->getright());
		return (root->value + tmp->value)/2;
	}

	return root->value;
}

/*==================================================================== AVLNode ======================================================================*/

template <typename T>
AVLNode<T>::AVLNode()
{
	value = 0;
	left = NULL;
	right = NULL;
	height = 0;
	count = 1;
}

template <typename T>
AVLNode<T>::AVLNode(T val)
{
	value = val;
	height = 0;
	count = 1;
	left = NULL;
	right = NULL;
}

template <typename T>
AVLNode<T>::AVLNode(T val, AVLNode<T> * left, AVLNode<T> * right)
{
	value = val;
	height = 0;
	count = 1;
	setleft(left);
	setright(right);
}

template <typename T>
AVLNode<T>::~AVLNode()
{
	value = 0;
	height = 0;
	count = 0;

	if (left != NULL)
	{
		delete left;
		left = NULL;
	}
	if (right != NULL)
	{
		delete right;
		right = NULL;
	}
}

template <typename T>
void AVLNode<T>::setleft(AVLNode<T> * node)
{
	left = node;
	int lhei = 0;
	int rhei = 0;
	int lcon = 0;
	int rcon = 0;
	if (left != NULL)
	{
		lhei = left->getheight()+1;
		lcon = left->getcount();
	}
	if (right != NULL)
	{
		rhei = right->getheight()+1;
		rcon = right->getcount();
	}
	height = lhei > rhei ? lhei : rhei;
	count = lcon + rcon + 1;
}

template <typename T>
void AVLNode<T>::setright(AVLNode<T> * node)
{
	right = node;
	int lhei = 0;
	int rhei = 0;
	int lcon = 0;
	int rcon = 0;
	if (left != NULL)
	{
		lhei = left->getheight()+1;
		lcon = left->getcount();
	}
	if (right != NULL)
	{
		rhei = right->getheight()+1;
		rcon = right->getcount();
	}
	height = lhei > rhei ? lhei : rhei;
	count = lcon + rcon + 1;
}

template class AVLTree<int>;
template class AVLTree<float>;
template class AVLTree<double>;