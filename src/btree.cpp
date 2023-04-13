/**
 * @author See Contributors.txt for code contributors and overview of BadgerDB.
 *
 * @section LICENSE
 * Copyright (c) 2012 Database Group, Computer Sciences Department, University of Wisconsin-Madison.
 */

#include "btree.h"
#include "filescan.h"
#include "exceptions/bad_index_info_exception.h"
#include "exceptions/bad_opcodes_exception.h"
#include "exceptions/bad_scanrange_exception.h"
#include "exceptions/no_such_key_found_exception.h"
#include "exceptions/scan_not_initialized_exception.h"
#include "exceptions/index_scan_completed_exception.h"
#include "exceptions/file_not_found_exception.h"
#include "exceptions/end_of_file_exception.h"
#include<vector>


// -----------------------------------------------------------------------------
// All the comments are documented in Integer Node handling functions
// since the algorithms for handling other 2 tpyes of keys are the same.
// -----------------------------------------------------------------------------

namespace badgerdb
{

// -----------------------------------------------------------------------------
// BTreeIndex::BTreeIndex -- Constructor
// -----------------------------------------------------------------------------

void BTreeIndex::initInt(std::string relationName){
	/* Psuedo code
	1. Initialize root node as NonLeafNode
	2. Initialize first two pages stored in root to ensure root node is always a NonLeafNode
	3. Initialize meta data 
	*/
	this->leafOccupancy = INTARRAYLEAFSIZE;
	this->nodeOccupancy = INTARRAYNONLEAFSIZE;
	PageId rootPageId;
	Page* rootPage;
	PageId metaPageId;
	Page* headerPage;
	Page* rp;
	PageId rpId;
	Page* lp;
	PageId lpId;
	this->bufMgr->allocPage(this->file, rootPageId, rootPage);
	this->bufMgr->allocPage(this->file, metaPageId, headerPage);
	this->bufMgr->allocPage(this->file, rpId, rp);
	this->bufMgr->allocPage(this->file, lpId, lp);

	// Initialize root node as NonLeafNode
	this->rootPageNum = rootPageId;
	NonLeafNodeInt *root = (NonLeafNodeInt*)rootPage;
	this->OriginalRootPageNum = rootPageNum;

	//Initialize first two pages stored in root to ensure root node is always a NonLeafNode
	root->size = 0;
	root->level = 1;
	root->pageNoArray[1] = rpId;
	root->pageNoArray[0] = lpId;
	((LeafNodeInt*)lp)->rightSibPageNo = rpId;

	//Initialize meta data 
	this->headerPageNum = metaPageId;
	IndexMetaInfo *metaInfo = (IndexMetaInfo*)headerPage;
	strcpy(metaInfo->relationName, relationName.c_str());
	metaInfo->attrByteOffset = this->attrByteOffset;
	metaInfo->attrType = this->attributeType;
	metaInfo->rootPageNo = this->rootPageNum;

	this->bufMgr->unPinPage(this->file, rootPageId, true);
	this->bufMgr->unPinPage(this->file, metaPageId, true);
	this->bufMgr->unPinPage(this->file, rpId, true);
	this->bufMgr->unPinPage(this->file, lpId, true);


}

void BTreeIndex::initDouble(std::string relationName){
	this->leafOccupancy = DOUBLEARRAYLEAFSIZE;
	this->nodeOccupancy = DOUBLEARRAYNONLEAFSIZE;
	PageId rootPageId;
	Page* rootPage;
	PageId metaPageId;
	Page* headerPage;
	Page* rp;
	PageId rpId;
	Page* lp;
	PageId lpId;
	this->bufMgr->allocPage(this->file, rootPageId, rootPage);
	this->bufMgr->allocPage(this->file, metaPageId, headerPage);
	this->bufMgr->allocPage(this->file, rpId, rp);
	this->bufMgr->allocPage(this->file, lpId, lp);

	this->rootPageNum = rootPageId;
	struct NonLeafNodeDouble *root = (NonLeafNodeDouble*)rootPage;
	this->OriginalRootPageNum = rootPageNum;
	root->size = 0;
	root->level = 1;
	
	root->pageNoArray[1] = rpId;
	root->pageNoArray[0] = lpId;
	((LeafNodeDouble*)lp)->rightSibPageNo = rpId;


	this->headerPageNum = metaPageId;
	IndexMetaInfo *metaInfo = (IndexMetaInfo*)headerPage;
	strcpy(metaInfo->relationName, relationName.c_str());
	metaInfo->attrByteOffset = this->attrByteOffset;
	metaInfo->attrType = this->attributeType;
	metaInfo->rootPageNo = this->rootPageNum;
	
	this->bufMgr->unPinPage(this->file, rootPageId, true);
	this->bufMgr->unPinPage(this->file, metaPageId, true);
	this->bufMgr->unPinPage(this->file, rpId, true);
	this->bufMgr->unPinPage(this->file, lpId, true);


}

void BTreeIndex::initString(std::string relationName){
	this->leafOccupancy = STRINGARRAYLEAFSIZE;
	this->nodeOccupancy = STRINGARRAYNONLEAFSIZE;
	PageId rootPageId;
	Page* rootPage;
	PageId metaPageId;
	Page* headerPage;
	Page* rp;
	PageId rpId;
	Page* lp;
	PageId lpId;
	this->bufMgr->allocPage(this->file, rootPageId, rootPage);
	this->bufMgr->allocPage(this->file, metaPageId, headerPage);
	this->bufMgr->allocPage(this->file, rpId, rp);
	this->bufMgr->allocPage(this->file, lpId, lp);

	this->rootPageNum = rootPageId;
	struct NonLeafNodeString *root = (NonLeafNodeString*)rootPage;
	this->OriginalRootPageNum = rootPageNum;
	root->size = 0;
	root->level = 1;

	
	root->pageNoArray[1] = rpId;
	root->pageNoArray[0] = lpId;
	((LeafNodeString*)lp)->rightSibPageNo = rpId;


	this->headerPageNum = metaPageId;
	IndexMetaInfo *metaInfo = (IndexMetaInfo*)headerPage;
	strcpy(metaInfo->relationName, relationName.c_str());
	metaInfo->attrByteOffset = this->attrByteOffset;
	metaInfo->attrType = this->attributeType;
	metaInfo->rootPageNo = this->rootPageNum;
	
	this->bufMgr->unPinPage(this->file, rootPageId, true);
	this->bufMgr->unPinPage(this->file, metaPageId, true);
	this->bufMgr->unPinPage(this->file, rpId, true);
	this->bufMgr->unPinPage(this->file, lpId, true);


}



BTreeIndex::BTreeIndex(const std::string & relationName,
		std::string & outIndexName,
		BufMgr *bufMgrIn,
		const int attrByteOffset,
		const Datatype attrType)
{
	/* Psuedo code
	1. Check if index file has already existed.
	if indexfile has already existed in current directory:
		open the existed one.
	else
		open a new one.
	2. Initialize the root node,the first two pages stored in root and the meta data
	3. Scan file to the end.
	*/

	std::ostringstream idxStr;
	idxStr << relationName << '.' << attrByteOffset;
	std::string indexName = idxStr.str();
	File* newFile;

	//Check if index file has already existed.
	bool indexAlreadyExisted = File::exists(indexName);
	if (indexAlreadyExisted)
		newFile = new BlobFile(indexName, false);
	else
		newFile = new BlobFile(indexName, true);
	this->file = newFile;
	this->attributeType = attrType;
	this->attrByteOffset = attrByteOffset;
	this->bufMgr = bufMgrIn;
	this->scanExecuting=false;
	outIndexName = indexName;
	switch(this->attributeType){
		case INTEGER:
		{
			//Initialize the root node,the first two pages stored in root and the meta data
			initInt(relationName);

			//Scan file to the end.
			FileScan *s = new FileScan(relationName, this->bufMgr);
			this->scanExecuting=true;
			try {
				RecordId rid;
				std::string record;
				void* key;
				while(1){
					s->scanNext(rid);
					record = s->getRecord();
					key = (void*)(record.c_str() + this->attrByteOffset);
					this->insertEntry(key, rid);
				}
			} catch (EndOfFileException e){
				delete(s);
				break;
			}
		}
		case DOUBLE:
		{
			initDouble(relationName);
			FileScan *s = new FileScan(relationName, this->bufMgr);
			this->scanExecuting=true;
			try {

				RecordId rid;
				std::string record;
				void* key;
				while(1){
					s->scanNext(rid);
					record = s->getRecord();
					key = (void*)(record.c_str() + this->attrByteOffset);

					this->insertEntry(key, rid);
				}
			} catch (EndOfFileException e){
				delete(s);
				break;
			}
		}
		case STRING:
		{

			initString(relationName);
			FileScan *s = new FileScan(relationName, this->bufMgr);
			this->scanExecuting=true;
			try {
				RecordId rid;
    			std::string record;
    			while(1){
				s->scanNext(rid);
				record = s->getRecord();
				std::string key = record.substr(this->attrByteOffset, 10);
				this->insertEntry(key.c_str(), rid);
				}
			} catch (EndOfFileException e){
				delete(s);
				break;
			}
		}
		default:
		{
			return;
		}
	}
}


// -----------------------------------------------------------------------------
// BTreeIndex::~BTreeIndex -- destructor
// -----------------------------------------------------------------------------

BTreeIndex::~BTreeIndex()
{
	this->bufMgr->cleanPages(this->file);
	bufMgr->flushFile(this->file);
	delete this->file;
}

// -----------------------------------------------------------------------------
// BTreeIndex::insertEntry
// -----------------------------------------------------------------------------

const void BTreeIndex::insertEntry(const void *key, const RecordId rid)
{
	switch (this-> attributeType)
	{
		case INTEGER:
		{
			//if root node is empty
			Page* curPage;
			bufMgr->readPage(file, this->rootPageNum, curPage);
			NonLeafNodeInt *rootNode = (NonLeafNodeInt *)curPage;
			if(rootNode->size == 0){
				int* tmp;
				tmp = (int*)key;
				rootNode->keyArray[0] = *tmp;
				rootNode->size += 1;
				rootNode ->level =1;
			}
			void * upkey;
			PageId upPid;
			insert(curPage, rootPageNum, false, key, rid, upkey , upPid);

			if(upkey != nullptr){
				Page* headerPage;
				Page* newPage;
				PageId newPid;

				bufMgr->allocPage(this->file,newPid,newPage);
				NonLeafNodeInt *newRoot = (NonLeafNodeInt *)newPage;
				int* tmp = (int *) upkey;
				newRoot->keyArray[0] = *(tmp);
				newRoot->pageNoArray[0] = this->rootPageNum;
				newRoot->pageNoArray[1] = upPid;
				newRoot->level = rootNode->level+1;
				newRoot->size = 1;
				bufMgr->unPinPage(file,rootPageNum,true);
				rootPageNum = newPid;
				bufMgr->readPage(file, headerPageNum, headerPage);
				((IndexMetaInfo*)headerPage)->rootPageNo = rootPageNum;
				bufMgr->unPinPage(file, rootPageNum,true);
			}
			return;

		}
		case DOUBLE:
		{
			//if root node is empty
			Page* curPage;
			bufMgr->readPage(file, this->rootPageNum, curPage);
			NonLeafNodeDouble *rootNode = (NonLeafNodeDouble *)curPage;
			if(rootNode->size == 0){
				double* tmp;
				tmp = (double*)key;
				rootNode->keyArray[0] = *tmp;
				rootNode->size += 1;
				rootNode ->level =1;
			}
			void * upkey;
			PageId upPid;
			insert(curPage, rootPageNum, false, key, rid, upkey , upPid);
			if(upkey != nullptr){
				Page* headerPage;
				Page* newPage;
				PageId newPid;

				bufMgr->allocPage(this->file,newPid,newPage);
				NonLeafNodeDouble *newRoot = (NonLeafNodeDouble *)newPage;
				double* tmp = (double *) upkey;
				newRoot->keyArray[0] = *(tmp);
				newRoot->pageNoArray[0] = this->rootPageNum;
				newRoot->pageNoArray[1] = upPid;
				newRoot->level = rootNode->level+1;
				newRoot->size = 1;
				bufMgr->unPinPage(file,rootPageNum,true);
				rootPageNum = newPid;
				bufMgr->readPage(file, headerPageNum, headerPage);
				((IndexMetaInfo*)headerPage)->rootPageNo = rootPageNum;
				bufMgr->unPinPage(file, rootPageNum,true);
			}
			return;

		}
		case STRING:
		{
			// std::cout<<"enter insert entry"<<std::endl;
			Page* curPage;
			// std::cout<<"root pid"<<rootPageNum<<std::endl;
			bufMgr->readPage(file, this->rootPageNum, curPage);
			// std::cout<<"enter insert entry 1."<<std::endl;
			NonLeafNodeString *rootNode = (NonLeafNodeString *)curPage;
			if(rootNode->size == 0){

				// std::cout<<"5"<<std::endl;
				strncpy(rootNode->keyArray[0] , (char*)key, 10);
				// std::cout<<*rootNode->keyArray[0]<<std::endl;
				rootNode->size += 1;
				rootNode->level = 1;
			}
			void * upkey;
			PageId upPid;
			// std::cout<<"3"<<std::endl;
			insert(curPage, rootPageNum, false, key, rid, upkey , upPid);
			// std::cout<<"4"<<std::endl;
			if(upkey != nullptr){
				Page* headerPage;
				Page* newPage;
				PageId newPid;
				bufMgr->allocPage(this->file,newPid,newPage);
				NonLeafNodeString *newRoot = (NonLeafNodeString *)newPage;

				strncpy(newRoot->keyArray[0], (char*)upkey, 10);
				newRoot->pageNoArray[0] = this->rootPageNum;
				newRoot->pageNoArray[1] = upPid;
				newRoot->level = rootNode->level+1;
				newRoot->size = 1;
				bufMgr->unPinPage(file,rootPageNum,true);
				rootPageNum = newPid;
				// std::cout<<"new root pid after split"<< newPid<<std::endl;
				bufMgr->readPage(file, headerPageNum, headerPage);
				((IndexMetaInfo*)headerPage)->rootPageNo = rootPageNum;
				bufMgr->unPinPage(file, rootPageNum,true);
			}

		}
		default:
		{
			return;
		}
	}
}

 void BTreeIndex::insert(Page *curPage, PageId curPageNum, bool isLeaf, const void *key, const RecordId rid, void *&upKey, PageId &upPid)
{
		/*
		1. If it is leaf node:
			1) determine if the node is full or not.
			if not full:
				insert the entry and split the node into two from the middle
				record the right sibling pid and key return back to the caller
			else:
				insert the entry into the correct position of node
		2. The node is non leaf node:
			1) find the the correct index (i) of the First Larger Index
			2) get page at PidArray[i] recursivelly calling the insert method until node level is 1
			3) if the level gets to 1, check the returned key and pid:
				if not null and node not full:
					insert the key and pid into the correct position of node
				else:
					insert the key and pid and split the node into two from the middle
					record the pid and key of the new node and return back to the caller
		*/
	switch(attributeType)
	{
		case INTEGER:
		{
			//reached leaf nodes
			LeafNodeInt *leaf = (LeafNodeInt *)curPage;
			if (isLeaf && leaf ->size < leafOccupancy) {
				//if not full insert
				insertIntLeaf(leaf,key,rid);
				upKey = nullptr;
				return;
			}
			if (isLeaf && leaf->size >= leafOccupancy){
				//if is full insert and split
				splitIntLeaf(leaf, key, rid, upKey, upPid);
				return;
			}
			//reached non-leaf nodes
			NonLeafNodeInt* currNode = (NonLeafNodeInt*)curPage;
			int firstLargerIndex = findTheFirstLargerIndexNonLeaf(curPage, key);
			PageId nextId = currNode->pageNoArray[firstLargerIndex];
			Page* nextPage;
			bufMgr->readPage(file, nextId, nextPage);
			bufMgr->unPinPage(file, curPageNum, false);
			//recusrively calling the method until level1
			if (currNode->level == 1){
				insert(nextPage,nextId, true, key, rid, upKey, upPid);
			}
			else{
				insert(nextPage,nextId, false, key, rid, upKey, upPid);
			}
			bufMgr->unPinPage(file, nextId, true);
			bufMgr->readPage(file, curPageNum, curPage);
			//check for the upkey and uppid 
			if(upKey != NULL && currNode->size < nodeOccupancy){
				//if not full, only insert
				insertIntNonLeaf(currNode, upKey, upPid);
				//bufMgr->unPinPage(file, currentPageNum, true);
				upKey = NULL;
				return;
			}
			if(upKey != NULL && currNode->size >= nodeOccupancy){
				void* tmpKey = upKey;
				PageId tmpupPid = upPid;
				//if is full insert and split
				splitIntNonLeaf(currNode, tmpKey, tmpupPid, upKey, upPid);
				return;

			}
			return;
		}
		case DOUBLE:
		{
			//reached leaf nodes
			LeafNodeDouble *leaf = (LeafNodeDouble *)curPage;
			if (isLeaf && leaf ->size < leafOccupancy) {
				insertDoubleLeaf(leaf,key,rid);
				upKey = nullptr;
				return;
			}
			if (isLeaf && leaf->size >= leafOccupancy){
				splitDoubleLeaf(leaf, key, rid, upKey, upPid);
				return;
			}
			//reached non-leaf nodes
			NonLeafNodeDouble* currNode = (NonLeafNodeDouble*)curPage;
			int firstLargerIndex = findTheFirstLargerIndexNonLeaf(curPage, key);
			PageId nextId = currNode->pageNoArray[firstLargerIndex];
			Page* nextPage;
			bufMgr->readPage(file, nextId, nextPage);
			bufMgr->unPinPage(file, curPageNum, false);
			if (currNode->level == 1)
				insert(nextPage,nextId, true, key, rid, upKey, upPid);
			else
				insert(nextPage,nextId, false, key, rid, upKey, upPid);
			bufMgr->unPinPage(file, nextId, true);

			bufMgr->readPage(file, curPageNum, curPage);
			if(upKey != NULL && currNode->size < nodeOccupancy){
				insertDoubleNonLeaf(currNode, upKey, upPid);
				//bufMgr->unPinPage(file, currentPageNum, true);
				upKey = NULL;
				return;
			}
			if(upKey != NULL && currNode->size >= nodeOccupancy){
				void* tmpKey = upKey;
				PageId tmpupPid = upPid;
				splitDoubleNonLeaf(currNode, tmpKey, tmpupPid, upKey, upPid);
				return;

			}
			return;
		}
		case STRING:
		{
			//reached leaf nodes
			LeafNodeString *leaf = (LeafNodeString *)curPage;
			if (isLeaf && leaf ->size < leafOccupancy) {
				// std::cout<<"22"<<std::endl;
				insertStringLeaf(leaf,key,rid);
				// std::cout<<"23"<<std::endl;
	// 				for(int i =  0; i < leaf->size; i++){
	// 	std::cout << leaf ->keyArray[i] <<"  " <<i <<"   ";
	// }
	// std::cout << std::endl;
				upKey = nullptr;
				return;
			}
			if (isLeaf && leaf->size >= leafOccupancy){
				// std::cout<<"24"<<std::endl;
				splitStringLeaf(leaf, key, rid, upKey, upPid);
				// std::cout<<"25"<<std::endl;
				return;
			}
			//reached non-leaf nodes
			NonLeafNodeString* currNode = (NonLeafNodeString*)curPage;
			int firstLargerIndex = findTheFirstLargerIndexNonLeaf(curPage, key);
			PageId nextId = currNode->pageNoArray[firstLargerIndex];
			Page* nextPage;
			bufMgr->readPage(file, nextId, nextPage);
			bufMgr->unPinPage(file, curPageNum, false);
			if (currNode->level == 1){
				// std::cout<<"21"<<std::endl;
				insert(nextPage,nextId, true, key, rid, upKey, upPid);
			}
			else{
				insert(nextPage,nextId, false, key, rid, upKey, upPid);
			}
			bufMgr->unPinPage(file, nextId, true);

			bufMgr->readPage(file, curPageNum, curPage);
			if(upKey != NULL && currNode->size < nodeOccupancy){
				// std::cout<<"28"<<std::endl;
				insertStringNonLeaf(currNode, upKey, upPid);
				// std::cout<<"29"<<std::endl;
				//bufMgr->unPinPage(file, currentPageNum, true);
				upKey = NULL;
				return;
			}
			if(upKey != NULL && currNode->size >= nodeOccupancy){
				void* tmpKey = upKey;
				PageId tmpupPid = upPid;
				splitStringNonLeaf(currNode, tmpKey, tmpupPid, upKey, upPid);
				return;

			}
			return;
		}
		default:
	
		break;
	}


}


void BTreeIndex::insertIntNonLeaf(NonLeafNodeInt * node, const void *key, const PageId pid){
	/*
	1) find the right place to insert
	2) move all larger item to the right by one position
	*/
	int i =  0;
	int tmpkey;
	PageId tmpid;
	for(i =  0; i < node->size; i++){
		//insert the key to coorect place, record the original value
		if(node ->keyArray[i] > *(int*) key){
			tmpkey = node ->keyArray[i];
			tmpid = node ->pageNoArray[i+1];
			node->keyArray[i] = *(int*) key;
			node->pageNoArray[i+1] = pid;
			break;
		}
	}

	// if kay is larger than all the keys in the node, then add to the end
	if(i == node->size){
		node->keyArray[node->size] = *(int*)key;
		node->pageNoArray[node->size + 1] = pid;
		node->size ++;
		return;
	}
	i++;

	// move larger item to the right
	for(; i< node-> size; i++){
		int tmpkey2 = node ->keyArray[i];
		PageId pid2 = node -> pageNoArray[i+1];

		node->keyArray[i] = tmpkey;
		node->pageNoArray[i+1] = tmpid;
		tmpkey = tmpkey2;
		tmpid = pid2;
	}

	node->keyArray[node->size] = tmpkey;
	node->pageNoArray[node->size +1] = tmpid;
	node->size++;
}

void BTreeIndex::insertDoubleNonLeaf(NonLeafNodeDouble * node, const void *key, const PageId pid){
	int i =  0;
	double tmpkey;
	PageId tmpid;
	for(i =  0; i < node->size; i++){
		//move larger key to the right
		if(node ->keyArray[i] > *(double*) key){
			tmpkey = node ->keyArray[i];
			tmpid = node ->pageNoArray[i+1];
			node->keyArray[i] = *(double*) key;
			node->pageNoArray[i+1] = pid;
			break;
		}
	}

	if(i == node->size){
		node->keyArray[node->size] = *(double*)key;
		node->pageNoArray[node->size + 1] = pid;
		node->size ++;
		return;
	}
	i++;

	for(; i< node-> size; i++){
		double tmpkey2 = node ->keyArray[i];
		PageId pid2 = node -> pageNoArray[i+1];

		node->keyArray[i] = tmpkey;
		node->pageNoArray[i+1] = tmpid;
		tmpkey = tmpkey2;
		tmpid = pid2;
	}
	node->keyArray[node->size] = tmpkey;
	node->pageNoArray[node->size +1] = tmpid;
	node->size++;
}

void BTreeIndex::insertStringNonLeaf(NonLeafNodeString * node, const void *key, const PageId pid){
	int i =  0;
	char* tmpkey;
	PageId tmpid;
	// std::cout<<"41"<<std::endl;
	for(i =  0; i < node->size; i++){
		//move larger key to the right
		if(strncmp(node ->keyArray[i] , (char*)key, 10)>0){
			// strncpy(tmpkey,node ->keyArray[i], 10);
			// tmpid = node ->pageNoArray[i+1];
			// strncpy(node->keyArray[i], (char*)key, 10);
			// node->pageNoArray[i+1] = pid;
			break;
		}
	}
	// std::cout<<"42"<<std::endl;
	for(int j = node ->size; j> i; j--){
		strncpy(node ->keyArray[j],node->keyArray[j-1],10);
		node->pageNoArray[j+1] = node->pageNoArray[j];
	}

	strncpy(node->keyArray[i],(char*)key,10);
	node->pageNoArray[i+1] = pid;
	// std::cout<<"13"<<std::endl;
	node->size++;

}

void BTreeIndex::insertIntLeaf(LeafNodeInt * leaf, const void *key, const RecordId rid){
	int i =  0;
	int tmpkey;
	RecordId tmpid;
	for(i =  0; i < leaf->size; i++){
		//move larger key to the right
		if(leaf ->keyArray[i] > *(int*) key){
			tmpkey = leaf ->keyArray[i];
			tmpid = leaf->ridArray[i];
			leaf->keyArray[i] = *(int*) key;
			leaf->ridArray[i] = rid;
			break;
		}
	}
	if(i == leaf-> size){
		leaf->keyArray[leaf->size] = *(int*)key;
		leaf->ridArray[leaf->size] = rid;
		leaf->size ++;
		return;
	}
	i++;
	for(; i< leaf-> size; i++){
		int tmpkey2 = leaf ->keyArray[i];
		RecordId rid2 = leaf -> ridArray[i];

		leaf->keyArray[i] = tmpkey;
		leaf->ridArray[i] = tmpid;
		tmpkey = tmpkey2;
		tmpid = rid2;
	}
	leaf->keyArray[leaf->size] = tmpkey;
	leaf->ridArray[leaf->size] = tmpid;
	leaf->size++;
}

void BTreeIndex::insertDoubleLeaf(LeafNodeDouble * leaf, const void *key, const RecordId rid){
	int i =  0;
	double tmpkey;
	RecordId tmpid;
	for(i =  0; i < leaf->size; i++){
		//move larger key to the right
		if(leaf ->keyArray[i] > *(double*) key){
			tmpkey = leaf ->keyArray[i];
			tmpid = leaf->ridArray[i];
			leaf->keyArray[i] = *(double*) key;
			leaf->ridArray[i] = rid;
			break;
		}
	}
	if(i == leaf-> size){
		leaf->keyArray[leaf->size] = *(double*)key;
		leaf->ridArray[leaf->size] = rid;
		leaf->size ++;
		return;
	}
	i++;
	for(; i< leaf-> size; i++){
		double tmpkey2 = leaf ->keyArray[i];
		RecordId rid2 = leaf -> ridArray[i];

		leaf->keyArray[i] = tmpkey;
		leaf->ridArray[i] = tmpid;
		tmpkey = tmpkey2;
		tmpid = rid2;
	}
	leaf->keyArray[leaf->size] = tmpkey;
	leaf->ridArray[leaf->size] = tmpid;
	leaf->size++;
}

void BTreeIndex::insertStringLeaf(LeafNodeString * leaf, const void *key, const RecordId rid){
	int i =  0;
	// std::cout<<"10"<<std::endl;
	RecordId tmpid;

	for(i =  0; i < leaf->size; i++){
		//move larger key to the right
		if(strncmp(leaf ->keyArray[i], (char*)key, 10)>0){
			break;
		}
	}
	for(int j = leaf ->size; j> i; j--){
		strncpy(leaf ->keyArray[j],leaf->keyArray[j-1],10);
		leaf->ridArray[j] = leaf->ridArray[j-1];
	}

	strncpy(leaf->keyArray[i],(char*)key,10);
	leaf->ridArray[i] = rid;
	// std::cout<<"13"<<std::endl;
	leaf->size++;
	// std::cout<<"15"<<std::endl;
}




void BTreeIndex::splitIntNonLeaf(NonLeafNodeInt* node, const void* key, PageId pid,  void *&upKey, PageId &upPid){
	/*
	1) record the smaller items and find the right place to insert
	2) append the insert entry
	3) append all the remaining larger items
	4) split from the middle, smaller half stores to the original node
	5) larger half stores into the new node
	*/
	Page* node2Page;
	PageId node2Id;
	this->bufMgr->allocPage(this->file, node2Id, node2Page);
	NonLeafNodeInt* node2 = (NonLeafNodeInt*) node2Page;
	node2 ->size = 0;
	int mid = (node->size +1) / 2;
	int newlen = node->size +1;
	int newKey[newlen];
	int newPid[newlen];
	int i;
	for(i =0; i< node->size; i++){
		if(*(int*)key < node->keyArray[i]){
			// find the right place to insert and record the entry
			newKey[i] = *(int*)key;
			newPid[i] = pid;
			break;
		}else{
			//record the smaller items
			newKey[i] = node->keyArray[i];
			newPid[i] = node->pageNoArray[i+1];
		}
	}
	// if key is larger than all the keys in the node, then add to the end
	if(i == node->size ){
		newKey[i] = *(int*) key;
		newPid[i] = pid;
	}else{
		//append all the remaining larger items
		for(; i< node->size ; i++){
			newKey[i+1] = node->keyArray[i];
			newPid[i+1] = node->pageNoArray[i+1];
		}
	}

	//smaller half stores to the original node
	for(int j = 0, k = 1; k<= mid ; k++, j++){
		node ->keyArray[j] = newKey[j];
		node ->pageNoArray[k] = newPid[j];
	}

	//larger half stores into the new node
	for(int j=0, k = mid+1; k < newlen; j++, k ++) {
		node2 ->keyArray[j] = newKey[k];
		node2 ->pageNoArray[j] = newPid[k-1];
	}
	node2 ->pageNoArray[newlen-mid-1] = newPid[newlen-1];
	node2->size = newlen - mid;
	node->size = mid;
	node2 ->level = node ->level;
	upKey = (void*) &newKey[mid];
	upPid = node2Id;
	this->bufMgr->unPinPage(file,node2Id,true);
}

void BTreeIndex::splitDoubleNonLeaf(NonLeafNodeDouble* node, const void* key, PageId pid,  void *&upKey, PageId &upPid){
	Page* node2Page;
	PageId node2Id;
	this->bufMgr->allocPage(this->file, node2Id, node2Page);
	NonLeafNodeInt* node2 = (NonLeafNodeInt*) node2Page;
	node2 ->size = 0;
	int mid = (node->size +1) / 2;
	int newlen = node->size +1;
	double newKey[newlen];
	int newPid[newlen];
	int i;
	for(i =0; i< node->size; i++){
		if(*(double*)key < node->keyArray[i]){
			newKey[i] = *(double*)key;
			newPid[i] = pid;
			break;
		}else{
			newKey[i] = node->keyArray[i];
			newPid[i] = node->pageNoArray[i+1];
		}
	}
	if(i == node->size ){
		newKey[i] = *(double*) key;
		newPid[i] = pid;
	}else{
		for(; i< node->size ; i++){
			newKey[i+1] = node->keyArray[i];
			newPid[i+1] = node->pageNoArray[i+1];
		}
	}

	for(int j = 0, k = 1; k<= mid ; k++, j++){
		node ->keyArray[j] = newKey[j];
		node ->pageNoArray[k] = newPid[j];
	}

	for(int j=0, k = mid+1; k < newlen; j++, k ++) {
		node2 ->keyArray[j] = newKey[k];
		node2 ->pageNoArray[j] = newPid[k-1];
	}
	node2 ->pageNoArray[newlen-mid-1] = newPid[newlen-1];
	node2->size = newlen - mid;
	node->size = mid;
	node2 ->level = node ->level;
	upKey = (void*) &newKey[mid];
	upPid = node2Id;
	this->bufMgr->unPinPage(file,node2Id,true);
}


void BTreeIndex::splitStringNonLeaf(NonLeafNodeString* node, const void* key, PageId pid,  void *&upKey, PageId &upPid){
	// std::cout<<"spliting" << std::endl; 
	Page* node2Page;
	PageId node2Id;
	this->bufMgr->allocPage(this->file, node2Id, node2Page);
	NonLeafNodeString* node2 = (NonLeafNodeString*) node2Page;
	node2 ->size = 0;
	int mid = (node->size +1) / 2;
	int newlen = node->size +1;
	char newKey[newlen][10];

	int newPid[newlen];
	int i;
	for(i =0; i< node->size; i++){
		if(strncmp((char*) key,node->keyArray[i], 10)<0){
			strncpy(newKey[i], (char*) key, 10);
			newPid[i] = pid;
			break;
		}else{
			strncpy(newKey[i],node->keyArray[i], 10);
			newPid[i] = node->pageNoArray[i+1];
		}
	}
	if(i == node->size ){
		strncpy(newKey[i],(char*) key, 120);
		newPid[i] = pid;
	}else{
		for(; i< node->size ; i++){
			strncpy(newKey[i+1],node->keyArray[i], 10);
			newPid[i+1] = node->pageNoArray[i+1];
		}
	}

	for(int j = 0, k = 1; k<= mid ; k++, j++){
		strncpy(node ->keyArray[j],newKey[j], 10);
		node ->pageNoArray[k] = newPid[j];
	}

	for(int j=0, k = mid+1; k < newlen; j++, k ++) {
		strncpy(node2 ->keyArray[j],newKey[k], 10);
		node2 ->pageNoArray[j] = newPid[k-1];
	}
	node2 ->pageNoArray[newlen-mid-1] = newPid[newlen-1];
	node2->size = newlen - mid;
	node->size = mid;
	node2 ->level = node ->level;
	// std::cout<<"splitStringNonLeaf"<<std::endl;

	const char* tmp = newKey[mid];
	upKey = (void*) tmp;
	// std::cout<<"splitStringNonLeaf finished."<<std::endl;
	upPid = node2Id;
	this->bufMgr->unPinPage(file,node2Id,true);

}


void BTreeIndex::splitIntLeaf(LeafNodeInt* node, const void* key, RecordId rid, void *&upKey, PageId &upPid){
	Page* node2Page;
	PageId node2Id;
	//std::cout<<"starting spliting leaf."<<std::endl;
	this->bufMgr->allocPage(this->file, node2Id, node2Page);
	LeafNodeInt* node2 = (LeafNodeInt*) node2Page;
	node2 ->size = 0;
	int mid = (node->size +1) / 2;
	int newlen = node->size +1;
	int newKey[newlen];
	RecordId newRid[newlen];
	int i;
	for(i =0; i< node->size; i++){
		if(*(int*)key < node->keyArray[i]){
			newKey[i] = *(int*)key;
			newRid[i] = rid;
			break;
		}else{
			newKey[i] = node->keyArray[i];
			newRid[i] = node->ridArray[i];
		}
	}
	if(i == node->size ){
		newKey[i] = *(int*) key;
		newRid[i] = rid;
	}else{
		for(; i< node->size ; i++){
			newKey[i+1] = node->keyArray[i];
			newRid[i+1] = node->ridArray[i];
		}
	}

	for(i = 0; i<mid; i++){
		node ->keyArray[i] = newKey[i];
		node ->ridArray[i] = newRid[i];
	}

//std::cout<<"before inserting node 2 "<<std::endl;
	for( int j=0, k = mid; k < newlen; j++, k ++) {
		node2 ->keyArray[j] = newKey[k];
		node2 ->ridArray[j] = newRid[k];
	}
	node2->size = newlen - mid;
	node->size = mid;
	//set right sibling
	node2->rightSibPageNo = node->rightSibPageNo;
	node->rightSibPageNo = node2Id;
	upKey = (void*) &(node2 ->keyArray[0]);
	upPid = node2Id;
		this->bufMgr->unPinPage(file,node2Id,true);
}

void BTreeIndex::splitDoubleLeaf(LeafNodeDouble* node, const void* key, RecordId rid, void *&upKey, PageId &upPid){
	Page* node2Page;
	PageId node2Id;
	this->bufMgr->allocPage(this->file, node2Id, node2Page);
	LeafNodeDouble* node2 = (LeafNodeDouble*) node2Page;
	node2 ->size = 0;
	int mid = (node->size +1) / 2;
	int newlen = node->size +1;
	double newKey[newlen];
	RecordId newRid[newlen];
	int i;
	for(i =0; i< node->size; i++){
		if(*(double*)key < node->keyArray[i]){
			newKey[i] = *(double*)key;
			newRid[i] = rid;
			break;
		}else{
			newKey[i] = node->keyArray[i];
			newRid[i] = node->ridArray[i];
		}
	}
	if(i == node->size ){
		newKey[i] = *(double*) key;
		newRid[i] = rid;
	}else{
		for(; i< node->size ; i++){
			newKey[i+1] = node->keyArray[i];
			newRid[i+1] = node->ridArray[i];
		}
	}

	for(i = 0; i<mid; i++){
		node ->keyArray[i] = newKey[i];
		node ->ridArray[i] = newRid[i];
	}

	for( int j=0, k = mid; k < newlen; j++, k++) {
		node2 ->keyArray[j] = newKey[k];
		node2 ->ridArray[j] = newRid[k];
	}
	node2->size = newlen - mid;
	node->size = mid;
	//set right sibling
	node2->rightSibPageNo = node->rightSibPageNo;
	node->rightSibPageNo = node2Id;
	upKey = (void*) &(node2 ->keyArray[0]);
	upPid = node2Id;
	this->bufMgr->unPinPage(file,node2Id,true);

}

void BTreeIndex::splitStringLeaf(LeafNodeString* node, const void* key, RecordId rid, void *&upKey, PageId &upPid){
	Page* node2Page;
	PageId node2Id;
	this->bufMgr->allocPage(this->file, node2Id, node2Page);
	LeafNodeString* node2 = (LeafNodeString*) node2Page;
	node2 ->size = 0;
	int mid = (node->size +1) / 2;
	int newlen = node->size +1;
	// std::cout<<"31"<<std::endl;
	char newKey[newlen][10];
	// std::cout<<"32"<<std::endl;
	RecordId newRid[newlen];
	int i;
	for(i =0; i< node->size; i++){
		if(strncmp((char*) key, node->keyArray[i], 10) < 0){
			// std::cout<<"33"<<std::endl;
			strncpy(newKey[i],(char*) key, 10);
			// std::cout<<"34"<<std::endl;
			newRid[i] = rid;
			break;
		}else{
			// std::cout<<"35"<<std::endl;
			strncpy(newKey[i],node->keyArray[i], 10);
			// std::cout<<"36"<<std::endl;
			newRid[i] = node->ridArray[i];
		}
	}
	if(i == node->size ){
		strncpy(newKey[i],(char*) key, 10);
		newRid[i] = rid;
	}else{
		for(; i< node->size ; i++){
			strncpy(newKey[i+1],node->keyArray[i], 10);
			newRid[i+1] = node->ridArray[i];
		}
	}

	for(i = 0; i<mid; i++){
		strncpy(node ->keyArray[i] , newKey[i], 10);
		node ->ridArray[i] = newRid[i];
	}

	for( int j=0, k = mid; k < newlen; j++, k ++) {
		strncpy(node2 ->keyArray[j] , newKey[k], 10);
		node2 ->ridArray[j] = newRid[k];
	}
	node2->size = newlen - mid;
	node->size = mid;
	//set right sibling
	node2->rightSibPageNo = node->rightSibPageNo;
	node->rightSibPageNo = node2Id;
	// std::cout<<"splitStringLeaf"<<std::endl;
	const char* tmp;
	tmp=node2 ->keyArray[0];
	upKey = (void*) tmp;
	// std::cout<<"splitStringLeaf finished."<<std::endl;
	upPid = node2Id;
	this->bufMgr->unPinPage(file,node2Id,true);

}





// -----------------------------------------------------------------------------
// BTreeIndex::startScan
// -----------------------------------------------------------------------------

const void BTreeIndex::startScan(const void* lowValParm,
				   const Operator lowOpParm,
				   const void* highValParm,
				   const Operator highOpParm)
{
	/* Psuedo Code


	1. Check if scan range is valid
	2. Conduct m-way search on B+ Tree
		(1) Traverse from root to the buttom NonLeafNode(level=1)
		(2) Find the leaf where the first larger leaf exists
	3. Check if the first larger key exists in keyArrays of leaf nodes
	if the first larger key not exists in current pointed node and
		current pointed node doesn't has the right sibling)
		Reach the very right end of the B+Tree
		--> lower bound of scanning range is larger the biggest key of the whole B+Tree
		--> no such key found
		return
	if (the first larger key not exists in current pointed node but
		current pointed node has the right sibling):
		move the current node pointer to the right sibling and 
		the key value stored in index 0 of keyArray should be the first larger key
		return
	return
	*/

	//Check if scan range is valid
	if((lowOpParm != GT && lowOpParm != GTE) || (highOpParm != LT && highOpParm != LTE)) {
		throw BadOpcodesException();
	}
	switch(this->attributeType){
		case INTEGER:
		{
			int lowerBound = *((int*)lowValParm);
			int higherBound = *((int*)highValParm);
			if(lowerBound > higherBound) throw BadScanrangeException();
			this->lowValInt = lowerBound;
			this->highValInt = higherBound;
			this->lowOp = lowOpParm;
			this->highOp = highOpParm;

			//Conduct m-way search on B+ Tree
			//Traverse from root to the buttom NonLeafNode(level=1)
			PageId curPageId = this->rootPageNum;
			Page* curPage;
			this->bufMgr->readPage(this->file, curPageId, curPage);
			NonLeafNodeInt* curPageNode = (NonLeafNodeInt*)curPage;
			while(curPageNode->level != 1) {
				int* key = &lowValInt;
				int firstLargerIndex = findTheFirstLargerIndexNonLeaf(curPage, (void*)key);
				struct NonLeafNodeInt* tmp = (NonLeafNodeInt*)(curPage);
				this->bufMgr->unPinPage(this->file, curPageId, false);
				PageId nextPageId = tmp->pageNoArray[firstLargerIndex];
				curPageId = nextPageId;
				this->bufMgr->readPage(this->file, curPageId, curPage);
				curPageNode = (NonLeafNodeInt*)curPage;
			}

			//Find the leaf where the first larger leaf exists
			int* key = &lowValInt;
			int firstLargerIndex = findTheFirstLargerIndexNonLeaf(curPage, (void*)key);
			struct NonLeafNodeInt* tmp = (NonLeafNodeInt*)(curPage);
			this->bufMgr->unPinPage(this->file, curPageId, false);

			//Check if the first larger key exists in keyArrays of leaf nodes
			PageId leafId = tmp->pageNoArray[firstLargerIndex];
			Page* leafPage;
			this->bufMgr->readPage(this->file, leafId, leafPage);
			LeafNodeInt* leafNode = (LeafNodeInt*)leafPage;
			bool found = false;
			for(int i=0;i<leafNode->size;i++){
				if (this->isLarger((void*)leafNode, i)){
					this->nextEntry = i;
					this->currentPageNum = leafId;
					this->currentPageData = leafPage;
					this->scanExecuting = true;
					found = true;
					break;
				}
			}
			/*
			if the first larger key not exists in current pointed node and
				current pointed node doesn't has the right sibling)
				Reach the very right end of the B+Tree
				--> lower bound of scanning range is larger the biggest key of the whole B+Tree
				--> no such key found
			*/
			if(!found &&  leafNode->rightSibPageNo== 0){
				this->scanExecuting = false;
				throw NoSuchKeyFoundException();
			}
			/*
			if (the first larger key not exists in current pointed node but
				current pointed node has the right sibling):
				move the current node pointer to the right sibling and 
				the key value stored in index 0 of keyArray should be the first larger key
			*/
			if(!found && leafNode->rightSibPageNo >= 0){
				Page* rightPage;
				this->bufMgr->readPage(this->file, leafNode->rightSibPageNo, rightPage);
				this->nextEntry = 0;
				this->currentPageNum = leafNode->rightSibPageNo;
				this->currentPageData = rightPage;
				this->scanExecuting = true;
				this->bufMgr->unPinPage(this->file, leafId, false);
			}
			//the first larger key found
			return;

		}
		case DOUBLE:
		{
			double lowerBound = *((double*)lowValParm);
			double higherBound = *((double*)highValParm);
			if(lowerBound > higherBound) throw BadScanrangeException();
			this->lowValDouble = lowerBound;
			this->highValDouble = higherBound;
			this->lowOp = lowOpParm;
			this->highOp = highOpParm;

			PageId curPageId = this->rootPageNum;
			Page* curPage;
			this->bufMgr->readPage(this->file, curPageId, curPage);
			NonLeafNodeDouble* curPageNode = (NonLeafNodeDouble*)curPage;

			while(curPageNode->level != 1) {
				double* key = &(this-> lowValDouble);
				int firstLargerIndex = findTheFirstLargerIndexNonLeaf(curPage, (void*)key);
				struct NonLeafNodeDouble* tmp = (NonLeafNodeDouble*)(curPage);
				this->bufMgr->unPinPage(this->file, curPageId, false);
				PageId nextPageId = tmp->pageNoArray[firstLargerIndex];
				curPageId = nextPageId;
				this->bufMgr->readPage(this->file, curPageId, curPage);
				curPageNode = (NonLeafNodeDouble*)curPage;
			}
			double* key = &(this-> lowValDouble);
			int firstLargerIndex = findTheFirstLargerIndexNonLeaf(curPage, (void*)key);
			struct NonLeafNodeDouble* tmp = (NonLeafNodeDouble*)(curPage);
			this->bufMgr->unPinPage(this->file, curPageId, false);
			PageId leafId = tmp->pageNoArray[firstLargerIndex];
			Page* leafPage;
			this->bufMgr->readPage(this->file, leafId, leafPage);
			LeafNodeDouble* leafNode = (LeafNodeDouble*)leafPage;

			bool found = false;
			for(int i=0;i<leafNode->size;i++){
				if (this->isLarger((void*)leafNode, i)){
					this->nextEntry = i;
					this->currentPageNum = leafId;
					this->currentPageData = leafPage;
					this->scanExecuting = true;
					found = true;
					break;
				}
			}
			if(!found && leafNode->rightSibPageNo == 0){
				this->scanExecuting = false;
				throw NoSuchKeyFoundException();
			}
			if(!found && leafNode->rightSibPageNo >= 0){
				Page* rightPage;
				this->bufMgr->readPage(this->file, leafNode->rightSibPageNo, rightPage);
				this->nextEntry = 0;
				this->currentPageNum = leafNode->rightSibPageNo;
				this->currentPageData = rightPage;
				this->scanExecuting = true;
				this->bufMgr->unPinPage(this->file, leafId, false);
			}
			return;
		}
		case STRING:
		{
			char* lowerBound = (char*)lowValParm;
			char* higherBound = (char*)highValParm;
			if( strncmp(lowerBound, higherBound, 10) > 0 ) throw BadScanrangeException();
			std::string lowStr = lowerBound;
			std::string highStr = higherBound;
			this->lowValString = lowStr.substr(0, 10);
			this->highValString = highStr.substr(0, 10);
			this->lowOp = lowOpParm;
			this->highOp = highOpParm;

			PageId curPageId = this->rootPageNum;
			Page* curPage;
			this->bufMgr->readPage(this->file, curPageId, curPage);
			NonLeafNodeString* curPageNode = (NonLeafNodeString*)curPage;

			while(curPageNode->level != 1) {
				const char* key = lowValString.c_str();
				int firstLargerIndex = findTheFirstLargerIndexNonLeaf(curPage, (void*)key);
				struct NonLeafNodeString* tmp = (NonLeafNodeString*)(curPage);
				this->bufMgr->unPinPage(this->file, curPageId, false);
				PageId nextPageId = tmp->pageNoArray[firstLargerIndex];
				curPageId = nextPageId;
				this->bufMgr->readPage(this->file, curPageId, curPage);
				curPageNode = (NonLeafNodeString*)curPage;
			}
			const char* key = lowValString.c_str();

			int firstLargerIndex = findTheFirstLargerIndexNonLeaf(curPage, (void*)key);
			struct NonLeafNodeString* tmp = (NonLeafNodeString*)(curPage);
			this->bufMgr->unPinPage(this->file, curPageId, false);
			PageId leafId = tmp->pageNoArray[firstLargerIndex];
			Page* leafPage;
			this->bufMgr->readPage(this->file, leafId, leafPage);
			LeafNodeString* leafNode = (LeafNodeString*)leafPage;

			bool found = false;
			for(int i=0;i<leafNode->size;i++){
				if (this->isLarger((void*)leafNode, i)){
					this->nextEntry = i;
					this->currentPageNum = leafId;
					this->currentPageData = leafPage;
					this->scanExecuting = true;
					found = true;
					break;
				}
			}
			if(!found && leafNode->rightSibPageNo == 0){
				this->scanExecuting = false;
				throw NoSuchKeyFoundException();
			}
			if(!found && leafNode->rightSibPageNo >= 0){
				Page* rightPage;
				this->bufMgr->readPage(this->file, leafNode->rightSibPageNo, rightPage);
				this->nextEntry = 0;
				this->currentPageNum = leafNode->rightSibPageNo;
				this->currentPageData = rightPage;
				this->scanExecuting = true;
				this->bufMgr->unPinPage(this->file, leafId, false);
			}
			return;
		}
		default:
		{

		}
	}

}

const bool BTreeIndex::isLarger(void *cur_node, int i){
	switch(this->attributeType){
		case INTEGER:
		{
			LeafNodeInt* leaf = ((LeafNodeInt*)cur_node);
			if(this->lowOp == GT)
				return leaf->keyArray[i] > this-> lowValInt;
			else
				return leaf->keyArray[i] >= this-> lowValInt;
		}
		case DOUBLE:
		{
			LeafNodeDouble* leaf = ((LeafNodeDouble*)cur_node);
			if(this->lowOp == GT)
				return leaf->keyArray[i] > this-> lowValDouble;
			else
				return leaf->keyArray[i] >= this-> lowValDouble;
		}
		case STRING:
		{
			LeafNodeString* leaf = ((LeafNodeString*)cur_node);
			if(this->lowOp == GT)
				return strncmp(leaf->keyArray[i] , (this-> lowValString).c_str(), 10) > 0;
			else
				return strncmp(leaf->keyArray[i] , (this-> lowValString).c_str(), 10) >= 0;
		}
		default:
		{
			return false;
		}
	}
	return false;
}


int BTreeIndex::findTheFirstLargerIndexNonLeaf(Page* curPage,const void* key)
{
	//Psuedo code: binary search 
	switch(this->attributeType){
		case INTEGER:
		{

			struct NonLeafNodeInt* tmp = (NonLeafNodeInt*)(curPage);
			int* k = (int*)key;
			int l = 0;
			int r = tmp->size-1;
			while(l<=r){
				int mid = (l+r)/2;

				if ( (tmp->keyArray[mid]) > *k)
					r = mid - 1;
				else
					l = mid + 1;
			}

			return l;
		}
		case DOUBLE:
		{
			struct NonLeafNodeDouble* tmp = (NonLeafNodeDouble*)(curPage);
			double* k = (double*)key;
			int l = 0;
			int r = tmp->size-1;
			while(l<=r){
				int mid = (l+r)/2;
				if ( (tmp->keyArray[mid]) > *(k))
					r = mid - 1;
				else
					l = mid + 1;
			}
			return l;
		}
		case STRING:
		{
			struct NonLeafNodeString* tmp = (NonLeafNodeString*)(curPage);
			int l = 0;
			int r = tmp->size-1;
			while(l<=r){
				int mid = (l+r)/2;
				if ( strncmp((tmp->keyArray[mid]), (char*)key, 10) > 0)
					r = mid - 1;
				else
					l = mid + 1;
			}
			return l;
		}
		default:
		{
			return INT_MAX;
		}
	}
}


// -----------------------------------------------------------------------------
// BTreeIndex::scanNext
// -----------------------------------------------------------------------------

const void BTreeIndex::scanNext(RecordId& outRid)
{
	/* Psuedo Code
	This API is only called after scan starts, meaning that both scan range and 
	lower bound are valid, so our goal is to use this API to iterate starting from 
	lower bound to higher bound of scanning range and find out "valid keys" in leaf nodes.
	[Note] def of "valid keys": keys that are not larger than higher bound and the largest 
	key of this B+Tree.

	1. Check if it is scanning
		if not --> throw exception
	2. Check if the nextEntry is out of bound(reach the right most end of B+Tree)
		if yes --> throw exception
	3. Check if cur_node->keyArray[nextEntry] < higher bound 
	if cur_node->keyArray[nextEntry] < higher bound:
		find the key!
		(1) store the key
		(2) nextEntry += 1(move nextEntry)
	4. Check if nextEntry has been moved:
	if nextEntry has been moved:
		if (nextEntry has reached the last position of keyArray and right sibling of 
		current exists):
			move the cur_node to its right sibling and set nextEntry as 0
	else:
		can't find more key --> scan complete (throw exception)
	return
	*/

	//Check if it is scanning
	if(!this->scanExecuting) throw ScanNotInitializedException();

	switch(this->attributeType){
		case INTEGER:
		{
			//Check if the nextEntry is out of bound(reach the right most end of B+Tree)
			LeafNodeInt *cur_node = (LeafNodeInt*)this->currentPageData;
			if (this->nextEntry >= cur_node -> size && cur_node->rightSibPageNo == 0) {
					throw IndexScanCompletedException();
			}

			//Check if cur_node->keyArray[nextEntry] < higher bound
			/*
			if cur_node->keyArray[nextEntry] < higher bound:
				find the key!
				(1) store the key
				(2) nextEntry += 1(move nextEntry)
			*/ 
			bool moved = false;
			if(this -> isSmaller((void*)cur_node)){
				moved = true;
				outRid = cur_node->ridArray[this->nextEntry];
				this->nextEntry ++;
			}

			//Check if nextEntry has been moved
			/*
			if nextEntry has been moved:
				if (nextEntry has reached the last position of keyArray and right sibling of 
				current exists):
					move the cur_node to its right sibling and set nextEntry as 0
			else:
				can't find more key --> scan complete (throw exception)
			*/
			if (moved){
				if (cur_node->rightSibPageNo > 0 && this->nextEntry >= cur_node -> size ){
					PageId tmp = this->currentPageNum;
					this->currentPageNum = cur_node->rightSibPageNo;
					this->bufMgr->unPinPage(this->file, tmp, false);
					this->bufMgr->readPage(this->file, this->currentPageNum, this->currentPageData);
					this->nextEntry = 0;
					cur_node = (LeafNodeInt*)this->currentPageData;
				}
				return;
			}else{
				throw IndexScanCompletedException();
			}
		}
		case DOUBLE:
		{
			LeafNodeDouble *cur_node = (LeafNodeDouble*)this->currentPageData;
			if (this->nextEntry >= cur_node -> size && cur_node->rightSibPageNo == 0) {
					throw IndexScanCompletedException();
			}
			bool moved = false;
			if(this -> isSmaller((void*)cur_node)){
				moved = true;
				outRid = cur_node->ridArray[this->nextEntry];
				this->nextEntry ++;
			}
			if (moved){
				if (this->nextEntry >= cur_node -> size && cur_node->rightSibPageNo > 0){
					PageId tmp = this->currentPageNum;
					this->currentPageNum = cur_node->rightSibPageNo;
					this->bufMgr->unPinPage(this->file, tmp, false);
					this->bufMgr->readPage(this->file, this->currentPageNum, this->currentPageData);
					this->nextEntry = 0;
					cur_node = (LeafNodeDouble*)this->currentPageData;
				}
				return;
			}else{
				throw IndexScanCompletedException();
			}
		}
		case STRING:
		{
			LeafNodeString *cur_node = (LeafNodeString*)this->currentPageData;
			if (this->nextEntry >= cur_node -> size && cur_node->rightSibPageNo == 0) {
					throw IndexScanCompletedException();
			}
			bool moved = false;
			if(this -> isSmaller((void*)cur_node)){
				moved = true;
				outRid = cur_node->ridArray[this->nextEntry];
				this->nextEntry ++;
			}
			if (moved){
				if (this->nextEntry >= cur_node -> size && cur_node->rightSibPageNo > 0){
					PageId tmp = this->currentPageNum;
					this->currentPageNum = cur_node->rightSibPageNo;
					this->bufMgr->unPinPage(this->file, tmp, false);
					this->bufMgr->readPage(this->file, this->currentPageNum, this->currentPageData);
					this->nextEntry = 0;
					cur_node = (LeafNodeString*)this->currentPageData;
				}
				return;
			}else{
				throw IndexScanCompletedException();
			}
		}
		default:
		{
			return;
		}
	}
}

const bool BTreeIndex::isSmaller(void *cur_node){
	switch(this->attributeType){
		case INTEGER:
		{
			LeafNodeInt* leaf = ((LeafNodeInt*)cur_node);
			if(this->highOp == LT)
				return leaf->keyArray[this->nextEntry] < this->highValInt;
			else
				return leaf->keyArray[this->nextEntry] <= this->highValInt;
		}
		case DOUBLE:
		{
			LeafNodeDouble* leaf = ((LeafNodeDouble*)cur_node);
			if(this->highOp == LT)
				return leaf->keyArray[this->nextEntry] < this->highValDouble;
			else
				return leaf->keyArray[this->nextEntry] <= this->highValDouble;
		}
		case STRING:
		{
			LeafNodeString* leaf = ((LeafNodeString*)cur_node);
			if(this->highOp == LT)
				return strncmp(leaf->keyArray[this->nextEntry] , (this->highValString).c_str(), 10) < 0;
			else
				return strncmp(leaf->keyArray[this->nextEntry], (this->highValString).c_str(), 10) <= 0;
		}
		default:
		{
			return false;
		}
	}
	return false;
}

// -----------------------------------------------------------------------------
// BTreeIndex::endScan
// -----------------------------------------------------------------------------
//
const void BTreeIndex::endScan() {
	if(!this->scanExecuting) throw ScanNotInitializedException();
	this->scanExecuting = false;
	this->bufMgr->unPinPage(this->file, this->currentPageNum, false);
}

}
