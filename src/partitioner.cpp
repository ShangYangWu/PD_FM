#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cassert>
#include <vector>
#include <cmath>
#include <map>
#include <climits>       // added 
#include "cell.h"
#include "net.h"
#include "partitioner.h"

using namespace std;


void Partitioner::parseInput(fstream& inFile)
{
    string str;
    // Set balance factor
    inFile >> str;
    _bFactor = stod(str);

    // Set up whole circuit
    while (inFile >> str) {
        if (str == "NET") {
            string netName, cellName, tmpCellName = "";
            inFile >> netName;
            int netId = _netNum;
            // cout << netName << " " << _netNum << endl;
            _netArray.push_back(new Net(netName));
            _netName2Id[netName] = netId;
            while (inFile >> cellName) {
                if (cellName == ";") {
                    tmpCellName = "";
                    break;
                }
                else {
                    // a newly seen cell
                    if (_cellName2Id.count(cellName) == 0) {
                        int cellId = _cellNum;
                        _cellArray.push_back(new Cell(cellName, 0, cellId));
                        _cellName2Id[cellName] = cellId;
                        _cellArray[cellId]->addNet(netId);
                        _cellArray[cellId]->incPinNum();
                        _netArray[netId]->addCell(cellId);
                        ++_cellNum;
                        tmpCellName = cellName;
                    }
                    // an existed cell
                    else {
                        if (cellName != tmpCellName) {
                            assert(_cellName2Id.count(cellName) == 1);
                            int cellId = _cellName2Id[cellName];
                            _cellArray[cellId]->addNet(netId);
                            _cellArray[cellId]->incPinNum();
                            _netArray[netId]->addCell(cellId);
                            tmpCellName = cellName;
                        }
                    }
                }
            }
            ++_netNum;
        }
    }
    return;
}


void Partitioner::printSummary() const
{
    cout << endl;
    cout << "==================== Summary ====================" << endl;
    cout << " Cutsize: " << _cutSize << endl;
    cout << " Total cell number: " << _cellNum << endl;
    cout << " Total net number:  " << _netNum << endl;
    cout << " Cell Number of partition A: " << _partSize[0] << endl;
    cout << " Cell Number of partition B: " << _partSize[1] << endl;
    cout << "=================================================" << endl;
    cout << endl;
    return;
}

void Partitioner::reportNet() const
{
    cout << "Number of nets: " << _netNum << endl;
    for (size_t i = 0, end_i = _netArray.size(); i < end_i; ++i) {
        cout << setw(8) << _netArray[i]->getName() << ": ";
        vector<int> cellList = _netArray[i]->getCellList();
        for (size_t j = 0, end_j = cellList.size(); j < end_j; ++j) {
            cout << setw(8) << _cellArray[cellList[j]]->getName() << " ";
        }
        cout << endl;
    }
    return;
}

void Partitioner::reportCell() const
{
    cout << "Number of cells: " << _cellNum << endl;
    for (size_t i = 0, end_i = _cellArray.size(); i < end_i; ++i) {
        cout << setw(8) << _cellArray[i]->getName() << ": ";
        vector<int> netList = _cellArray[i]->getNetList();
        for (size_t j = 0, end_j = netList.size(); j < end_j; ++j) {
            cout << setw(8) << _netArray[netList[j]]->getName() << " ";
        }
        cout << endl;
    }
    return;
}

void Partitioner::writeResult(fstream& outFile)
{
    stringstream buff;
    buff << _cutSize;
    outFile << "Cutsize = " << buff.str() << '\n';
    buff.str("");
    buff << _partSize[0];
    outFile << "G1 " << buff.str() << '\n';
    for (size_t i = 0, end = _cellArray.size(); i < end; ++i) {
        if (_cellArray[i]->getPart() == 0) {
            outFile << _cellArray[i]->getName() << " ";
        }
    }
    outFile << ";\n";
    buff.str("");
    buff << _partSize[1];
    outFile << "G2 " << buff.str() << '\n';
    for (size_t i = 0, end = _cellArray.size(); i < end; ++i) {
        if (_cellArray[i]->getPart() == 1) {
            outFile << _cellArray[i]->getName() << " ";
        }
    }
    outFile << ";\n";
    return;
}

void Partitioner::clear()
{
    for (size_t i = 0, end = _cellArray.size(); i < end; ++i) {
        delete _cellArray[i];
    }
    for (size_t i = 0, end = _netArray.size(); i < end; ++i) {
        delete _netArray[i];
    }
    return;
}

void Partitioner::partition()
{
    /*set balance bond*/
    _bond = (1-_bFactor)/2*_cellNum;;

    /* init partition */
    initPart(0);

    /* iterate */
    while(1){
        // init tracing record
        _moveStack.clear();

        // iterate
        ++_iterNum;
        while(_unlockNum[0]+_unlockNum[1]!=0){
            iterate();
        }

        // for all partial sum of _maxAccGain=0 then stop iterating
        if(_maxAccGain == 0){
            // recover to the best
            for(int round=_moveNum-1; round>_bestMoveNum-1; --round){ 
                const int cellId = _moveStack[round];
                const bool part = _cellArray[cellId]->getPart();
                moveCell(cellId, part);
                // cout << "trace back to moving "<< _cellArray[cellId]->getName() << endl;
            }
            _cutSize -= _maxAccGain;

            cout << "Total iterations: "<< _iterNum << endl;
            // printSummary();
            return;
        }

        // reset all cells
        for(const auto &item : _cellArray){
            item->unlock();
            ++_unlockNum[item->getPart()];
            item->setGain(0);
        }

        // recover to the best
        for(int round=_moveNum-1; round>_bestMoveNum-1; --round){ 
            const int cellId = _moveStack[round];
            const bool part = _cellArray[cellId]->getPart();
            moveCell(cellId, part);
            // cout << "trace back to moving "<< _cellArray[cellId]->getName() << endl;
        }
        _cutSize -= _maxAccGain;
        // cout << _iterNum << endl;
        // printSummary();

        initPart(1);
    }
}

void Partitioner::initPart(const bool iter){
    if(iter == 0){
        // step1: place cells in each party balanced
        int balance = _cellNum/2;
        for(size_t i=0; i < balance; ++i){
            ++_partSize[1];
            ++_unlockNum[1];
            _cellArray[i]->setPart(1);
            for(const auto &item : _cellArray[i]->getNetList()){
                _netArray[item]->incPartCount(1);
            }
        }
        for(size_t i=balance; i < _cellNum; ++i){
            ++_partSize[0];
            ++_unlockNum[0];
            // _cellArray[i]->setPart(0); // already initiated
            for(const auto &item : _cellArray[i]->getNetList()){
                _netArray[item]->incPartCount(0);
            }
        }

        // init cutSize
        for(const auto &item : _netArray){
            if(item->getPartCount(0)>0 && item->getPartCount(1)>0){
                ++_cutSize;
            }
        }

        // init iterNum
        _iterNum = 0;
    }

    // step2: initiate gain
    for(const auto &item : _netArray){
        for(const auto &it : item->getCellList()){           
            Cell* const cell = _cellArray[it];
            const bool from = cell->getPart();
            const int FromCount = item->getPartCount(from);
            const int ToCount = item->getPartCount(!from);

            // From = 1 => Gain++
            if(FromCount == 1){
                cell->incGain();
            }
            // To = 0 => Gain--
            if(ToCount == 0){
                cell->decGain();
            }  
        }

    }

    // build bList
    for(const auto &item : _cellArray){
        Node* const it = item->getNode();
        const int gain = item->getGain();

        if(_maxGainCell == NULL){
            _maxGainCell = it;
        }
        addNode(it, item->getPart(), gain);
        if(gain >= _cellArray[_maxGainCell->getId()]->getGain()){
            _maxGainCell = it;
        }
    }

    // init iter para
    _moveNum = 0;
    _bestMoveNum = 0;
    _maxAccGain = INT_MIN;

    // printBList();
    // printSummary();
}

void Partitioner::iterate(){
    // step1: decide FromSide and lock _maxGainCell
    const int maxGainCellId = _maxGainCell->getId();
    Cell* const maxGainCell = _cellArray[maxGainCellId];
    const int maxGainCellGain = _cellArray[maxGainCellId]->getGain();
    const bool From = _cellArray[maxGainCellId]->getPart();

    // lock maxGainCell
    rmNode(_maxGainCell, From, maxGainCellGain);
    maxGainCell->lock();
    // cout << _cellArray[_maxGainCell->getId()]->getName() << " is locked: " << maxGainCell->getLock() << endl;
    --_unlockNum[From];

    // update _accGain and _maxAccGain in this iteration
    if(_moveNum == 0){
        _accGain =  maxGainCellGain;
    }else{
        _accGain +=  maxGainCellGain;
    }
    ++_moveNum;
    if(_accGain >= _maxAccGain){
        _maxAccGain = _accGain;
        _bestMoveNum = _moveNum;
    }

    // step2: update gain for each node before move
    // before _maxGainCell moves => ToCount=0:gain++ / ToCount=1:gain(to)--
    for(const auto &item : maxGainCell->getNetList()){
        const int FromCount = _netArray[item]->getPartCount(From);
        const int ToCount = _netArray[item]->getPartCount(!From);
        if(ToCount == 0){
            // update bList[all] && ToCount=0:gain++
            for(const auto &it : _netArray[item]->getCellList()){
                Cell* const cell = _cellArray[it];
                if(cell->getLock() == 0){
                    Node* const node = cell->getNode();
                    const bool party = cell->getPart();
                    rmNode(node, party, cell->getGain());
                    cell->incGain();
                    addNode(node, party, cell->getGain());
                }
            }
        }
        if(ToCount == 1){
            // update bList[!From] && ToCount=1:gain(to)--
            for(const auto &it : _netArray[item]->getCellList()){
                Cell* const cell = _cellArray[it];
                if(cell->getPart() == (!From) && cell->getLock() == 0){
                    Node* const node = cell->getNode();
                    const bool party = cell->getPart();
                    rmNode(node, party, cell->getGain());
                    cell->decGain();
                    addNode(node, party, cell->getGain());
                }
            }
        }
    }

    // step3: move _maxGainCell
    moveCell(maxGainCellId, From);
    _moveStack.emplace_back(maxGainCellId);

    // step4: update gain for each node before move
    // After _maxGainCell moves => FromCount=0:gain-- / FromCount=1:gain(from)++
    for(const auto &item : maxGainCell->getNetList()){
        const int FromCount = _netArray[item]->getPartCount(From);
        const int ToCount = _netArray[item]->getPartCount(!From);
        if(FromCount == 0){
            // update bList[all] && FromCount=0:gain--
            for(const auto &it : _netArray[item]->getCellList()){
                Cell* const cell = _cellArray[it];
                if(cell->getLock() == 0){
                    Node* node = cell->getNode();
                    bool party = cell->getPart();
                    rmNode(node, party, cell->getGain());
                    cell->decGain();
                    addNode(node, party, cell->getGain());
                }
            }
        }
        if(FromCount == 1){
            // update bList[From] && FromCount=1:gain(from)++
            for(const auto &it : _netArray[item]->getCellList()){
                Cell* const cell = _cellArray[it];
                if(cell->getPart() == From && cell->getLock() == 0){
                    Node* node = cell->getNode();
                    bool party = cell->getPart();
                    rmNode(node, party, cell->getGain());
                    cell->incGain();
                    addNode(node, party, cell->getGain());
                }
            }
        }
    }

    // step5: select new _maxGainCell for the next iteration
    const int part0SizeAftermoved = _partSize[0]-1;
    const int part1SizeAftermoved = _partSize[1]-1;
    const bool bucket0Empty = _bList[0].empty();
    const bool bucket1Empty = _bList[1].empty();
    // no candidates
    if(bucket0Empty && bucket1Empty){
        // cout <<"G0 && G1 are both empty"<< endl;
        // printBList();
        // printSummary();
        return;
    }else{
        if(bucket1Empty){
            // check if G0 is movable
            if(part0SizeAftermoved >= _bond){
                // only G0 is movable
                map<int, Node*>::reverse_iterator maxBucketGain0 = _bList[0].rbegin();
                
                // maxGain in _bList[0] for maxGain
                _maxGainCell = maxBucketGain0->second;
                // cout <<"max in G0 "<< _cellArray[_maxGainCell->getId()]->getName() << " " << _cellArray[_maxGainCell->getId()]->getGain() << endl;
            }else{
                // G0 is unmovable
                // cout <<"G1 is empty, but G0 is unmovable"<< endl;
                // printBList();
                return;
            }

        }else{
            if(bucket0Empty){
                // check if G1 is movable
                if(part1SizeAftermoved >= _bond){
                    // only G1 is movable
                    map<int, Node*>::reverse_iterator maxBucketGain1 = _bList[1].rbegin();

                    // maxGain in _bList[1] for maxGain
                    _maxGainCell = maxBucketGain1->second;
                    // cout <<"max in G1 "<< _cellArray[_maxGainCell->getId()]->getName() << " " << _cellArray[_maxGainCell->getId()]->getGain() << endl;
                }else{
                    // cout <<"G0 is empty, but G1 is unmovable"<< endl;
                    // printBList();
                    // printSummary();
                    return;
                }

            }else{
                // 4 conditions
                if(part0SizeAftermoved >= _bond && part1SizeAftermoved >= _bond){
                    // G1 G2 are both movable
                    map<int, Node*>::reverse_iterator maxBucketGain0 = _bList[0].rbegin(); 
                    map<int, Node*>::reverse_iterator maxBucketGain1 = _bList[1].rbegin();

                    // compare and chose _maxGainCell
                    if(maxBucketGain0->first >= maxBucketGain1->first){
                        _maxGainCell = maxBucketGain0->second;
                    }
                    else{
                        _maxGainCell = maxBucketGain1->second;
                    }
                    // cout <<"max "<< _cellArray[_maxGainCell->getId()]->getName() << " " << _cellArray[_maxGainCell->getId()]->getGain() << endl;
                }else{
                    if(part0SizeAftermoved >= _bond){
                        // only G0 is movable
                        map<int, Node*>::reverse_iterator maxBucketGain0 = _bList[0].rbegin();
                        
                        // maxGain in _bList[0] for maxGain
                        _maxGainCell = maxBucketGain0->second;
                        // cout <<"max in G0 "<< _cellArray[_maxGainCell->getId()]->getName() << " " << _cellArray[_maxGainCell->getId()]->getGain() << endl;
                    }else{
                        if(part1SizeAftermoved >= _bond){
                            // only G1 is movable
                            map<int, Node*>::reverse_iterator maxBucketGain1 = _bList[1].rbegin();

                            // maxGain in _bList[1] for maxGain
                            _maxGainCell = maxBucketGain1->second;
                            // cout <<"max in G1 "<< _cellArray[_maxGainCell->getId()]->getName() << " " << _cellArray[_maxGainCell->getId()]->getGain() << endl;
                        }else{
                            // cout <<"G0 && G1 are not empty, but are both not movable"<< endl;
                            // printBList();
                            // printSummary();
                            return;
                        }
                    }
                }
            }
        }
    }

    // printBList();
    // printSummary();
}

void Partitioner::addNode(Node* const node, const bool party, const int gain){
    // bList -> [node]
    if(getBList(party)[gain] == NULL){
        _bList[party][gain] = node;
        node->setNext(NULL);
        node->setPrev(NULL);
        // cout << "add bList["<< party <<"]["<< gain <<"] -> [" << _cellArray[node->getId()]->getName() << "]"<< endl;
        return;
    }
    // bList -> [node] -> node
    else{
        node->setNext(_bList[party][gain]);
        node->getNext()->setPrev(node);
        node->setPrev(NULL);
        _bList[party][gain] = node;
        // cout << "add bList["<< party <<"]["<< gain <<"] -> [" << _cellArray[node->getId()]->getName() << "] -> " << _cellArray[node->getNext()->getId()]->getName() << endl;
        return;
    }
}

void Partitioner::rmNode(Node* const node, const bool party, const int gain){
    Node* const next = node->getNext();
    Node* const prev = node->getPrev();

    // bList -> node -> [node]
    if(next == NULL && prev != NULL){
        prev->setNext(NULL);
        node->setPrev(NULL);
        // cout << "remove bList["<< party <<"]["<< gain <<"] -> node -> [" << _cellArray[node->getId()]->getName() << "]" << endl;
        return;
    }    
    // bList -> node -> [node] -> node
    else if(next != NULL && prev != NULL){
        prev->setNext(next);
        next->setPrev(prev);
        node->setNext(NULL);
        node->setPrev(NULL);
        // cout << "remove bList["<< party <<"]["<< gain <<"] -> node -> [" << _cellArray[node->getId()]->getName() << "] -> node" << endl;
        return;
    }
    // bList -> [node] -> node
    else if(next != NULL ){
        _bList[party][gain] = next;
        next->setPrev(NULL);
        node->setNext(NULL);
        // cout << "remove bList["<< party <<"]["<< gain <<"] -> [" << _cellArray[node->getId()]->getName() << "] -> node" << endl;
        return;
    }   
    // bList -> [node]
    else{
		_bList[party].erase(gain);
        // cout << "remove bList["<< party <<"]["<< gain <<"] -> [" << _cellArray[node->getId()]->getName() << "]" << endl;
        return;
    }
}

void Partitioner::moveCell(const int id, const bool part){
    // update cutlist and net status
    Cell* const cell = _cellArray[id];
    for(const auto &item : cell->getNetList()){
        Net* const net = _netArray[item];
        net->decPartCount(part);
        net->incPartCount(!part);
    }
    // update partSize and move cell
    --_partSize[part];
    cell->move();
    ++_partSize[!part];
}

void Partitioner::printBList(){
    cout << "\n------------------------------" << endl;
    cout << "<" << " iterate " << _iterNum  << " - " << _moveNum << " >" << endl;
    cout << "Blist 0" << endl;
    cout << getPartSize(0) << endl;
    for(const auto &item : _bList[0]){
        Node* node = item.second;
        cout << item.first << " ";
        while(node!=NULL){
            cout << _cellArray[node->getId()]->getName() << " ";
            node = node->getNext();
        }
        cout << endl;
    }
    cout << "\nBlist 1" << endl;
    cout << getPartSize(1) << endl;
    for(const auto &item : _bList[1]){
        Node* node = item.second;
        cout << item.first << " ";
        while(node!=NULL){
            cout << _cellArray[node->getId()]->getName() << " ";
            node = node->getNext();
        }
        cout << endl;
    }
    cout << "\nMaxGainCell: "<< _cellArray[_maxGainCell->getId()]->getName() << " Gain:" << _cellArray[_maxGainCell->getId()]->getGain() << endl;
    cout << "AccGain: " << _accGain << " / " << "MaxGain: " << _maxAccGain << endl; 
    cout << "------------------------------" << endl;
}
