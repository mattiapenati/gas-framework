/*
 * Copyright (c) 2008, Davide Ferrarese & Mattia Penati
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. Neither the name of the <ORGANIZATION> nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef _GAS_LIST_H_
#define _GAS_LIST_H

#include <cstddef>

namespace Common {
	/* A node of the list */
	template<typename T>
	class Node {
		private:
			T val_;
			Node *next_;
			Node *prev_;
		public:
			Node(T val): val_(val), next_(0), prev_(0) {};
			void Next(Node *next) { next_ = next; };
			void Prev(Node *prev) { prev_ = prev; };
	};
	/* A list */
	template<typename T>
	class List {
		private:
			Node<T> *head_;
			Node<T> *current_;
			size_t L_;
			Node<T> *Move(size_t const);
		public:
			List();
			~List();
			inline T &operator()(size_t const);
			inline T const &operator()(size_t const) const;
			void Insert(T const, size_t const);
			void Remove(size_t const);
			inline size_t Size();
			inline bool Empty();
	};

	template<typename T>
	Node<T> *List<T>::Move(size_t const i) {
		assert(i < L_);
		Node<T> *t;
		t = head_;
		if (i < L_/2) {
			range(k, 0, i) t = t->next_;
		} else {
			range(k, 0, L_-i) t = t->prev_;
		}
		return t;
	};

	template<typename T>
	List<T>::List(): head_(0), current_(0), L_(0) {
	};

	template<typename T>
	T &List<T>::operator()(size_t const i) {
		return Move(i)->val_;
	};

	template<typename T>
	T const &List<T>::operator()(size_t const i) const {
		return Move(i)->val_;
	}

	template<typename T>
	void List<T>::Insert(T val, size_t const i) {
		assert(i <= L_);
		if (i == L_) {
			if (L_ == 0) {
				head_ = new Node<T>(val);
				head_->prev_ = head_;
				head_->next_ = prev_;
			} else {
				head_->prev_->next_ = new Node<T>(val);
				head_->prev_->next_->prev_ = head_->prev_;
				head_->prev_ = head_->prev_->next_;
				head_->prev_->next_ = head_;
			}
		} else {
			Node<T> *t;
			t = Move(i);
			t->prev_->next_ = new Node<T>(val);
			t->prev_->next_->prev_ = t->prev_;
			t->prev_ = t->prev_->next_;
			t->prev_->next_ = t;
			if (i == 0)
				head_ = t->prev_;
		}
		++L;
	}

	template<typename T>
	void List<T>::Remove(size_t const i) {
		Node<T> *t;
		t = Move(i);
		t->next_->prev_ = t->prev_;
		t->prev_->next_ = t->next_;
		if (i == 0)
			head_ = t->next_;
		delete t;
		--L_;
	}

	template<typename T>
	size_t List<T>::Size() {
		return L_;
	}

	template<typename T>
	bool List<T>::Empty() {
		return L_ == 0;
	}
}

#endif
