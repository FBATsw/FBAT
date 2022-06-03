// MySimplePointerList.h
#pragma once

template <class T> class MySimplePointerList {
	
	public:
	T*	node;
	MySimplePointerList *next;
	
	MySimplePointerList(T *p1=NULL, MySimplePointerList *p2=NULL) { node = p1; next = NULL; if (p2) { while (p2->next) p2=p2->next; p2->next = this; } }
	~MySimplePointerList() { if (node) delete node; if (next) delete next; }
	MySimplePointerList* detach(MySimplePointerList* ptr);
	MySimplePointerList* detach(T* pnode);
	void detachnode() { for (MySimplePointerList *lst=this; lst!=NULL; lst=lst->next) lst->node = NULL; }
};

template <class T> MySimplePointerList<T>* MySimplePointerList<T>::detach(MySimplePointerList<T>* ptr) {
	
	if (ptr==NULL)
		return this;
		
	MySimplePointerList* lst;
	
	if (ptr==this) {
		lst = next;
		next = NULL;
		return lst;
	}
	
	for (lst=this; lst->next && lst->next!=ptr; lst=lst->next)
		;
	if (lst->next==ptr) {
		lst->next = ptr->next;
		ptr->next = NULL;
	}
	
	return this;

}

template <class T> MySimplePointerList<T>* MySimplePointerList<T>::detach(T* pnode) {
	
	MySimplePointerList* lst, ptr;
	
	if (pnode==node) {
		lst = next;
		next = NULL;
		return lst;
	}
	
	for (lst=this; lst->next && lst->next->node!=pnode; lst=lst->next)
		;
	if (lst->next!=NULL) {
		ptr = lst->next;
		lst->next = ptr->next;
		ptr->next = NULL;
	}
	
	return this;

}
