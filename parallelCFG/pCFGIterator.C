#include "pCFGIterator.h"
#include <stack>
#include <vector>
#include <fstream>

#define HEX(x) setw(2) << std::hex << (int)x

int pCFGIteratorDebugLevel = 2;
vector<std::stack<DataflowNode> > pCFGNode::nodestack;

struct state{
    pCFGNode pcfg_node;
    set<unsigned int> blocked;
    set<unsigned int> active;
    set<unsigned int> released;

    state(){

    }

    state(pCFGNode p, set<unsigned int> b, set<unsigned int> a, set<unsigned int> r){
        this->pcfg_node = p;
        this->blocked = b;
        this->active = a;
        this->released = r;
    }
};

std::stack<state> globalstate;

void pCFGIterator::genInitState(const Function& func, const pCFGNode& n, const NodeState& state,
                      vector<Lattice*>& initLattices, vector<NodeFact*>& initFacts)
{                  
}

void pCFGIterator::copyPSetState(const Function& func, const pCFGNode& n, 
                   unsigned int srcPSet, unsigned int tgtPSet, NodeState& state,
                   vector<Lattice*>& lattices, vector<NodeFact*>& facts, 
                   ConstrGraph* partitionCond, bool omitRankSet)
{
    assert(0);
}

bool pCFGIterator::initPSetDFfromPartCond(const Function& func, const pCFGNode& n, unsigned int pSet,
                                const vector<Lattice*>& dfInfo, const vector<NodeFact*>& facts,
                                ConstrGraph* partitionCond)
{
    assert(0);
    return false;
}

void pCFGIterator::mergePCFGStates(const list<unsigned int>& pSetsToMerge, const pCFGNode& n, const Function& func,
                         NodeState& staet, const vector<Lattice*>& dfInfo, map<unsigned int, unsigned int>& pSetMigrations)
{
    assert(0);
}
    
void pCFGIterator::matchSendsRecvs(const pCFGNode& n, const vector<Lattice*>& dfInfo, NodeState* state, 
                     // Set by analysis to identify the process set that was split
                     unsigned int& splitPSet,
                     vector<ConstrGraph*>& splitConditions, 
                      vector<DataflowNode>& splitPSetNodes,
                     // for each split process set, true if its active and false if it is blocked
                     vector<bool>&         splitPSetActive,
                     // for each process set that was not split, true if becomes active as a result of the match,
                     // false if its status doesn't change
                     vector<bool>&         pSetActive,
                     const Function& func, NodeState* fState)
{
    assert(0);
}

// transfer function
// decides if the given pSet is blocked or split or dead
// if needs to split, pushes the descendants into splitPSetNodes
// if split -> pushes the ConstrGraph corresponding to that into splitConditions
bool pCFGIterator::transfer(const pCFGNode& n, 
                            unsigned int pSet, 
                            const Function& func,                  
                            NodeState& state, 
                            const vector<Lattice*>& dfInfo,
                            bool& deadPSet, 
                            bool& splitPSet, 
                            vector<DataflowNode>& splitPSetNodes,
                            bool& splitPNode, 
                            vector<ConstrGraph*>& splitConditions, 
                            bool& blockPSet)
{
    assert(0);
    return false;
}

// transfer function
// decides if the given pSet is blocked or split or dead
// if needs to split, pushes the descendants into splitPSetNodes
// if split -> pushes the ConstrGraph corresponding to that into splitConditions
bool pCFGIterator::transfer(const pCFGNode& n, 
                            unsigned int pSet, 
                            const Function& func,                  
                            NodeState& state, 
                            const vector<Lattice*>& dfInfo,
                            bool& isDeadPSet, 
                            bool& isSplitPSet, 
                            vector<DataflowNode>& splitPSetNodes,
                            bool& isSplitPNode,                             
                            bool& isBlockPSet,
                            bool& isMergePSet)
{
    bool modified = false;
    // Get the ROSE_VisitorPattern instance
    boost::shared_ptr<IntraPCFGTransferVisitor> 
        transferVisitor = boost::shared_ptr<IntraPCFGTransferVisitor> 
        (new pCFGIteratorTransfer (n, pSet, func, state, dfInfo, isDeadPSet, isSplitPSet, splitPSetNodes, isSplitPNode, isBlockPSet, isMergePSet, this->mda));

    // get the node on which visitor pattern needs to applied
    const DataflowNode& dfNode = n.getCurNode(pSet);

    // get the node
    SgNode* sgn = dfNode.getNode();

    // set the handler
    sgn->accept(*transferVisitor);

    modified = transferVisitor->finish() || modified;

    return modified;    
}

void pCFGIterator::resetPSet(unsigned int pSet, vector<Lattice*>& dfInfo)
{
    assert(0);
}


// pCFG_FWDataflow::runAnalysis invokes runAnalysis_pCFG
// Walks the pCFG to create dot representation
// Annotations used to match send/recv calls, process set names
// Transfer functions are only used to determine
//     1. NextState for a processSet
//     2. to block a particular processSet
//     3. Split a process set based on rank dependent condition
bool pCFGIterator::runAnalysis_pCFG(const Function& func, NodeState* state, pCFG_Checkpoint* chkpt)
{
	clock_t t;
    t = clock();
    std::cout << "CLIENT RUN ANALYSIS PCFG\n";
    string indent="";
    DataflowNode funcCFGStart = cfgUtils::getFuncStartCFG(func.get_definition());
    DataflowNode funcCFGEnd = cfgUtils::getFuncEndCFG(func.get_definition());

    ostringstream funcNameAsString;
    funcNameAsString << "Function " << func.get_name().getString() << "()";   
    if(pCFGIteratorDebugLevel >= 1) {
        Dbg::enterFunc(funcNameAsString.str());
        Dbg::dbg << indent << "Entering : " << funcNameAsString.str() << endl;
    }

    // current pCFG node
    pCFGNode curNode;

    // currently active
    set<unsigned int> activePSets;
    // currently blocked
    set<unsigned int> blockedPSets;
    // currently resumed from blocked
    // skip transfer function
    set<unsigned int> releasedPSets;

    // psets requiring merge
    set<unsigned int> mergePSets;

    // end sets
    set<unsigned int> endPSets;

    // Not restarting from a checkpoint
    // start with single process set
    // initialize to start from top of function
    if(chkpt == NULL) {
        vector<DataflowNode> initDFNodes;
        initDFNodes.push_back(funcCFGStart);
        curNode.init(initDFNodes);
        stack<DataflowNode> temp;
        vector<DataflowNode> temp2;
        (curNode.nodestack).push_back(temp);
        (curNode.sendSets).push_back(temp2);
        curNode.nodestack[0].push(funcCFGStart);
        activePSets.insert(0);

        // set up dot properties for initial node
        // curNode.set_id(1);
//        curNode.set_dot_attr("blue");
        // curNode.set_dot_attr(0);
        // nodeouts << curNode.toDot() << endl;
        npcfgnodes = 1;
    }
    // restarting from chkpt
    else {
        curNode = chkpt->n;
        activePSets = chkpt->activePSets;
        blockedPSets = chkpt->blockedPSets;
        releasedPSets = chkpt->releasedPSets;
        cout << indent << "RESTARTING FROM CHKPT " << chkpt->str(indent+"    ") << endl;
    }

    // bool shouldMatch = false;
    bool movedPSet = false;
    // outer-loop
    // apply dataflow on all process sets
    // match send/recv when all process sets are blocked
    ofstream fout("matches.txt");
do
{
    while(activePSets.size() > 0)
    {
        bool modified = false;
        // move until every pset is blocked
        // perform send/recv matching after this loop

        while(activePSets.size() > 0)
        {
            // cout << "activePsets.size() =  "<< activePSets.size() << " loop..\n";
            int curPSet;
            // The current process set is the next active one or the process set that 
            // performed the split from which we're restarting 
            // (the latter special condition is only needed to keep the output more readable)
            if(chkpt==NULL) {
                // pick in canonical order
                curPSet = *(activePSets.begin());
            }
            else {
                curPSet = chkpt->splitPSet;
            }

            // printSubGraphInit(curPSet);
                            
            // If we're restarting, we don't need the checkpoint any more
            if(chkpt) {
                delete chkpt;
                chkpt = NULL;
            }
            // process the curPset until it is blocked
            while(1) {
                if(curPSet >= curNode.nodestack.size()){
                    stack<DataflowNode> temp;
                    curNode.nodestack.push_back(temp);
                    vector<DataflowNode> temp2;
        			(curNode.sendSets).push_back(temp2);
                }
                NodeState* state = NULL;
                pCFGNode descNode(curNode);
                if(releasedPSets.erase(curPSet) > 0) {
                    descNode.advanceOut(curPSet);
                    if(movedPSet) {
                        movedPSet = false;                            
                    }
                    else {
                    }
                    modified = true;
                }
                else {
                    vector<Lattice*> dfInfoNewBelow;                                        
                    bool isSplitPSet = false;
                    bool isSplitPNode = false;
                    bool isBlockPSet = false;
                    bool isDeadPSet = false;
                    bool isMergePSet = false;
                   vector<DataflowNode> splitPSetNodes;
                    modified = transfer(descNode, curPSet, func, *state, dfInfoNewBelow,
                                        isDeadPSet, isSplitPSet, splitPSetNodes, isSplitPNode, isBlockPSet, isMergePSet);
                    // cout << "transfer end\n";

                    if(curNode.getCurNode(curPSet) == funcCFGEnd) {
                        // break
                        // move curPSet from active to end set
                        movePSet(curPSet, endPSets, activePSets);
                        break;
                    }

                    if(isDeadPSet) {
                        // Cannot be in this node
                        // stop progress
                        return false;
                    }
                    // if the analysis needs to split
                    // condition independent of 
                    else if(isSplitPNode) {

                    }
                    //NOTE:descNode should be updated by transfer function
                    else if(isSplitPSet) {
                        // cout << "isSplitPSet\n";
                        modified = true;
                        performPSetSplit(curPSet, curNode, descNode, splitPSetNodes, activePSets);
                        
                    }
                    // if curPSet wants to block
                    else if(isBlockPSet) {
                        // cout << "isBlockPSet\n";
                        movePSet(curPSet, blockedPSets, activePSets);
                        if(!movedPSet) {
                            // curNode.updateAncsNode(curPSet, curNode);
                        }
                        else {
                            movedPSet = false;
                        }
                        movedPSet = true;
                        printSubGraphEnd();
                        break;
                    }
                    else if(isMergePSet) {                        
                        // cout << "isMergePSet\n";
                        movePSet(curPSet, mergePSets, activePSets);
                        if(!movedPSet) {
                            // curNode.updateAncsNode(curPSet, curNode);
                        }
                        else {
                            movedPSet = false;
                        }
                        movedPSet = true;                       
                        // printSubGraphEnd();
                        break;
                    }
                    else {
                        // cout << "Dbg::dbg\n";
                        Dbg::dbg << indent << descNode.str() <<"\n";
                        descNode.advanceOut(curPSet);
                        // cout<<"2\n";                        
                        if(movedPSet) {
                            // do not update ancestor as it was set earlier for this pset
                            // reset the flag
                            movedPSet = false;                            
                        }
                        else {
                            // cout<<"not movedPSet\n";
                            // descNode.updateAncsNode(curPSet, curNode);
                        }
                        modified = true;
                    }
                    // <<<<<<<<<<< TRANSFER FUNCTION >>>>>>>>>>>>
                }// end else                               

                if(modified) {
                    curNode = descNode;
                    modified = false;
                }
                // this processSet is blocked
                //NOTE: modified = false -> fixpoint reached
                else {
                    return true;
                }                
                
            } // end of while(1) -> curPset is blocked or reached end of function
                        
            //movedPSet = true;
        } // end of while(activePsets.size > 0) -> all psets are blocked or require merge

        if(mergePSets.size() > 1) {
            unsigned int activepset = *(mergePSets.begin());
            DataflowNode curDfNode = curNode.getCurNode(activepset); // make a copy
            // vector<DataflowNode> initDF;
            // initDF.push_back(curDfNode);

            pCFGNode descNode(curNode.pSetDFNodes, curNode.visited, curNode.prev_send, curNode.prev_recv, curNode.matchings, curNode.sendSets);

            // cout << "\ndescNode->pSetDFNodes.size() = "<< descNode.pSetDFNodes.size() << endl;
            curNode = descNode;

            
            // mark visited
            visitedPCFGNodes.insert(descNode);
            
            releasedPSets.insert(activepset);
            activePSets.insert(activepset);
            mergePSets.clear();
            // we dont want to process blocked yet
        }

        else if(blockedPSets.size() > 0) {
            set<unsigned int> sendPSets;
            set<unsigned int> recvPSets;
            // cout << "blockPSets.size() = " << blockedPSets.size() << endl;
            bool check = filterSendRecv(curNode, blockedPSets, sendPSets, recvPSets, activePSets, releasedPSets);
            if(check){
		        ROSE_ASSERT(sendPSets.size() + recvPSets.size() == blockedPSets.size());

		        bool matched = matchSendRecv(curNode, sendPSets, recvPSets, blockedPSets, activePSets, releasedPSets);
		        if(matched) {
		            cout << "MATCH.....\n";
		        }
		        else{
		            fout << "\nPOTENTIAL DEADLOCK!!!\n";
		            fout.close();
		            cout << "Not matched!\n";
					cout << "\n time = " << (float)(clock() - t) / CLOCKS_PER_SEC << endl;
		            return true;
		        }
		    }
        }            
        //<<<<<< MatchSendRecv >>>>>>>>>

    }  // end of do-while
    
    //TODO: Check if all pSets are at funcCFGEnd
    //TODO: Create chkpts if not present
    //TODO: split(set<chkpts>)     // inserts new chkpts to set of already existing chkpts
    // runAnalysis() takes care of starting from the chkpt
    // All process sets are blocked and they cannot be unblocked via send-receive matching.
        
    // Check if this state was reached because all process sets are at the end of the application 
    int curPSet=0;
    vector<DataflowNode>::const_iterator it;
    for(it=curNode.getPSetDFNodes().begin(); 
        it!=curNode.getPSetDFNodes().end(); it++, curPSet++)
    {
        // If some process set is not at the end of this function
        if(*it != funcCFGEnd)
        {
            // Dbg::dbg << indent << "WARNING: in un-releaseable state process set "<<curPSet<<" is not at function end! Instead, n.pSetDFNodes["<<curPSet<<"]="<<(*it).str()<<endl;
            //ROSE_ASSERT(0);
            // Checkpoint this analysis. We may return to it later and discover that not all these process states
            // are actually possible.
            // if(analysisDebugLevel>=1) 
            //     Dbg::dbg << indent << "@@@ Shelving this blocked partition for now.\n";
            // set<pCFG_Checkpoint*> chkpt;
            // chkpt.insert(new pCFG_Checkpoint(curNode, func, state, activePSets, blockedPSets));
            // split(chkpt);
            // cout << "\nReached deadlock!!\n\n";
            break;
        }
    }
    if(it == curNode.getPSetDFNodes().end()){
        vector<match> all_matches = curNode.matchings;
        fout << "***************\n";
        for(int i=0;i<all_matches.size();i++){
            fout << "sender = (" << all_matches[i].sender << ", " << all_matches[i].send_id << ")\nreceiver = (" << all_matches[i].receiver << ", " << all_matches[i].recv_id << ")\n\n";
        }
        fout << "\n\n";
    }
    // cout<<"Program finished!!\n";
    if(globalstate.empty())
        break;
    struct state cool;
    cool = globalstate.top();
    globalstate.pop();
    // cout<<"Stack popped\n";
    curNode = cool.pcfg_node;
    blockedPSets = cool.blocked;
    activePSets = cool.active;
    releasedPSets = cool.released;

    set<unsigned int> sendPSets;
    set<unsigned int> recvPSets;
    // cout << "blockPSets.size() = " << blockedPSets.size() << endl;
    bool check = filterSendRecv(curNode, blockedPSets, sendPSets, recvPSets, activePSets, releasedPSets);
    if(check){
	    ROSE_ASSERT(sendPSets.size() + recvPSets.size() == blockedPSets.size());

	    bool matched = matchStar(curNode, sendPSets, recvPSets, blockedPSets, activePSets, releasedPSets);
    }

} while(!globalstate.empty());
    
        
    // At this point the analysis has reached the end of the function and does not need to propagate
    // any dataflow further down this function's pCFG. We will now return without creating a checkpoint.
    fout.close();    
    if(analysisDebugLevel>=1) Dbg::exitFunc(funcNameAsString.str());
    printSubGraphEnd();
	cout << "\n time = " << (float)(clock() - t) / CLOCKS_PER_SEC << endl;
        
    return true;
}

// update descNode
void pCFGIterator::performPSetSplit(int curPSet, const pCFGNode& curNode,
                                    pCFGNode& descNode, vector<DataflowNode>& splitPSetNodes, set<unsigned int>& activePSets)
{
    // remove curPSet from the pcfg node
    descNode.removePSet(curPSet);
    int _psetindex;
    stack<DataflowNode> temp = descNode.nodestack[curPSet];
    descNode.nodestack.erase(descNode.nodestack.begin() + curPSet);

    vector<DataflowNode> temp2 = descNode.sendSets[curPSet];
    descNode.sendSets.erase(descNode.sendSets.begin() + curPSet);
    
    vector<DataflowNode>::iterator dfI;
    for(dfI = splitPSetNodes.begin(); dfI != splitPSetNodes.end(); dfI++) {
        _psetindex = descNode.createPSet(*dfI);
        descNode.nodestack.push_back(temp);
        descNode.sendSets.push_back(temp2);
        // insert _psetindex into activeSets
        //TODO: need some checks here
        activePSets.insert(_psetindex);
        // descNode.updateAncsNode(_psetindex, curNode);
    }
}

bool pCFGIterator::filterSendRecv(pCFGNode& pcfgn, set<unsigned int>& blockedPSets, 
                                  set<unsigned int>& sendPSets,
                                  set<unsigned int>& recvPSets,
                                  set<unsigned int>& activePSets,
                                  set<unsigned int>& releasedPSets)
{
    set<unsigned int>::iterator it;
    for(it = blockedPSets.begin(); it != blockedPSets.end(); it++) {
        const DataflowNode& dfnode = pcfgn.getCurNode(*it);
        SgNode* sgn = dfnode.getNode();
        Function callee(isSgFunctionCallExp(sgn));

        if(callee.get_name().getString() == "MPI_Send") {
            sendPSets.insert(*it);
            DataflowNode dfsend = pcfgn.getCurNode(*it);
            ((pcfgn.sendSets)[*it]).push_back(dfsend);
            movePSet(*it, activePSets, blockedPSets);
	        releasedPSets.insert(*it);
	        return 0;
	        break;
        }
        else if(callee.get_name().getString() == "MPI_Recv") {
            recvPSets.insert(*it);
        }        
    }
    return 1;
}

// One match at a time
bool pCFGIterator::matchSendRecv(pCFGNode& pcfgn,
                                 set<unsigned int>& sendPSets,
                                 set<unsigned int>& recvPSets,
                                 set<unsigned int>& blockedPSets,
                                 set<unsigned int>& activePSets,
                                 set<unsigned int>& releasedPSets)
{
    set<unsigned int>::iterator send_it, recv_it;
    bool matched = false;
    // cout<<"entered matchSendRecv\n";

    // for(recv_it = recvPSets.begin(); recv_it != recvPSets.end(); recv_it++) {
        for(int r=0;r<recvPSets.size();r++){
        recv_it = recvPSets.begin();
        for(int i=0;i<r;i++){
            recv_it++;
        }

        // try to match each recv for a given send
        const DataflowNode& dfRecv = pcfgn.getCurNode(*recv_it);
        SgNode* sgnRecv = dfRecv.getNode();
        pcfgMatchAnnotation *recvAnnotation = dynamic_cast<pcfgMatchAnnotation*> (sgnRecv->getAttribute("pCFGAnnotation"));
        ROSE_ASSERT(recvAnnotation != NULL);
        pair<string, int> recvSource = recvAnnotation->getSource();
        pair<string, int> recvTarget = recvAnnotation->getTarget();
        // for(send_it = sendPSets.begin(); send_it != sendPSets.end(); send_it++) {
            for(int i=0;i<pcfgn.sendSets.size();i++){
            for(int j=0;j<pcfgn.sendSets[i].size();j++){
            const DataflowNode& dfSend = pcfgn.sendSets[i][j];
            SgNode* sgnSend = dfSend.getNode();
            pcfgMatchAnnotation *sendAnnotation = dynamic_cast<pcfgMatchAnnotation*> (sgnSend->getAttribute("pCFGAnnotation"));
            ROSE_ASSERT(recvAnnotation != NULL);
            pair<string, int> sendSource = sendAnnotation->getSource();
            pair<string, int> sendTarget = sendAnnotation->getTarget();

            if(recvSource.second <= 0){
                if(recvSource.first == sendTarget.first){
                    pCFGNode newpcfg = pcfgn;
                    newpcfg.prev_send = i;
                    newpcfg.prev_recv = r;
                    set<unsigned int> blocked = blockedPSets;
                    set<unsigned int> active = activePSets;
                    set<unsigned int> released = releasedPSets;

                    state st(newpcfg, blocked, active, released);
                    globalstate.push(st);
                    // cout<<"\nStack pushed (matchSendRecv)!!!\n\n";
                    matched = true;
                    // cout << "sendSource = " << sendSource.first << ", " << sendSource.second << endl;
                    // cout << "recvTarget = " << recvTarget.first << ", " << recvTarget.second << endl;
                    match m(sendSource.first, sendSource.second, recvTarget.first, recvTarget.second);
                    (pcfgn.matchings).push_back(m);
                    pcfgn.sendSets[i].erase(pcfgn.sendSets[i].begin() + j);
                    i = pcfgn.sendSets.size();
                    break;
                }
            }

            if(sendSource.first == recvTarget.first && recvSource.first == sendTarget.first) {
                pCFGNode newpcfg = pcfgn;
                newpcfg.prev_send = i;
                newpcfg.prev_recv = r;
                set<unsigned int> blocked = blockedPSets;
                set<unsigned int> active = activePSets;
                set<unsigned int> released = releasedPSets;

                state st(newpcfg, blocked, active, released);
                globalstate.push(st);
                // cout<<"\nStack pushed (matchSendRecv)!!!\n\n";
                
                // we have a match
                matched = true;
                // cout << "sendSource = " << sendSource.first << ", " << sendSource.second << endl;
                // cout << "recvSource = " << recvSource.first << ", " << recvSource.second << endl;
                match m(sendSource.first, sendSource.second, recvSource.first, recvSource.second);
                (pcfgn.matchings).push_back(m);
                pcfgn.sendSets[i].erase(pcfgn.sendSets[i].begin() + j);
                break;
            }            
            }
        }

        if(matched)
            break;
    }

    if(matched) {
        // int sendPSet = *send_it;
        int recvPSet = *recv_it;
        // move send/ recv psets from blocked to active
        // pCFGNode& ref = const_cast<pCFGNode&>(pcfgn); // downgrade the const cast
        // printEdge(ref.getCurAncsNode(sendPSet), ref.getCurAncsNode(recvPSet), "red");
        // movePSet(sendPSet, activePSets, blockedPSets);
        movePSet(recvPSet, activePSets, blockedPSets);
        // releasedPSets.insert(sendPSet);
        releasedPSets.insert(recvPSet);        
    }
    else {
        cout << "ERROR: No matches available in blocked state\n";
    }

    return matched;
}

bool pCFGIterator::matchStar(pCFGNode& pcfgn,
                                 set<unsigned int>& sendPSets,
                                 set<unsigned int>& recvPSets,
                                 set<unsigned int>& blockedPSets,
                                 set<unsigned int>& activePSets,
                                 set<unsigned int>& releasedPSets)
{
    set<unsigned int>::iterator send_it, recv_it;
    bool matched = false;
    // cout<<"entered matchStar\n";
    // for(recv_it = recvPSets.begin(); recv_it != recvPSets.end(); recv_it++) {
    // cout<<"recvPSets.size() = " << recvPSets.size() << endl;
    // cout<<"pcfgn.prev_recv = " << pcfgn.prev_recv << endl;
        for(int r=pcfgn.prev_recv;r<recvPSets.size();r++){
        recv_it = recvPSets.begin();
        for(int i=0;i<r;i++){
            recv_it++;
        }

        // try to match each recv for a given send
        const DataflowNode& dfRecv = pcfgn.getCurNode(*recv_it);
        SgNode* sgnRecv = dfRecv.getNode();
        pcfgMatchAnnotation *recvAnnotation = dynamic_cast<pcfgMatchAnnotation*> (sgnRecv->getAttribute("pCFGAnnotation"));
        ROSE_ASSERT(recvAnnotation != NULL);
        pair<string, int> recvSource = recvAnnotation->getSource();
        pair<string, int> recvTarget = recvAnnotation->getTarget();
        // for(send_it = sendPSets.begin(); send_it != sendPSets.end(); send_it++) {
            for(int i=pcfgn.prev_send+1;i<pcfgn.sendSets.size();i++){
            for(int j=0;j<pcfgn.sendSets[i].size();j++){
            const DataflowNode& dfSend = pcfgn.sendSets[i][j];
            SgNode* sgnSend = dfSend.getNode();
            pcfgMatchAnnotation *sendAnnotation = dynamic_cast<pcfgMatchAnnotation*> (sgnSend->getAttribute("pCFGAnnotation"));
            ROSE_ASSERT(recvAnnotation != NULL);
            pair<string, int> sendSource = sendAnnotation->getSource();
            pair<string, int> sendTarget = sendAnnotation->getTarget();

            if(recvSource.second <= 0){
                if(recvSource.first == sendTarget.first){
                    pCFGNode newpcfg = pcfgn;
                    newpcfg.prev_send = i;
                    newpcfg.prev_recv = r;
                    set<unsigned int> blocked = blockedPSets;
                    set<unsigned int> active = activePSets;
                    set<unsigned int> released = releasedPSets;

                    state st(newpcfg, blocked, active, released);
                    globalstate.push(st);
                    // cout<<"\nStack pushed (matchStar)!!!\n\n";
                    matched = true;
                    // cout << "sendSource = " << sendSource.first << ", " << sendSource.second << endl;
                    // cout << "recvTarget = " << recvTarget.first << ", " << recvTarget.second << endl;
                    match m(sendSource.first, sendSource.second, recvTarget.first, recvTarget.second);
                    (pcfgn.matchings).push_back(m);
                    pcfgn.sendSets[i].erase(pcfgn.sendSets[i].begin() + j);
                    break;
                }
            }

            if(sendSource.first == recvTarget.first && recvSource.first == sendTarget.first) {
                pCFGNode newpcfg = pcfgn;
                newpcfg.prev_send = i;
                newpcfg.prev_recv = r;
                set<unsigned int> blocked = blockedPSets;
                set<unsigned int> active = activePSets;
                set<unsigned int> released = releasedPSets;

                state st(newpcfg, blocked, active, released);
                globalstate.push(st);
                // cout<<"\nStack pushed (matchSendRecv)!!!\n\n";
                    
                // we have a match
                matched = true;
                // cout << "sendSource = " << sendSource.first << ", " << sendSource.second << endl;
                // cout << "recvSource = " << recvSource.first << ", " << recvSource.second << endl;
                match m(sendSource.first, sendSource.second, recvSource.first, recvSource.second);
                (pcfgn.matchings).push_back(m);
                pcfgn.sendSets[i].erase(pcfgn.sendSets[i].begin() + j);
                break;
            }            
            }
        }

        if(matched)
            break;
        pcfgn.prev_send=-1;
    }

    if(matched) {
        // int sendPSet = *send_it;
        int recvPSet = *recv_it;
        // move send/ recv psets from blocked to active
        // pCFGNode& ref = const_cast<pCFGNode&>(pcfgn); // downgrade the const cast
        // printEdge(ref.getCurAncsNode(sendPSet), ref.getCurAncsNode(recvPSet), "red");
        // movePSet(sendPSet, activePSets, blockedPSets);
        movePSet(recvPSet, activePSets, blockedPSets);
        // releasedPSets.insert(sendPSet);
        releasedPSets.insert(recvPSet);        
    }
    else {
        cout << "ERROR: No matches available in blocked state\n";
    }

    return matched;
}


void pCFGIterator::printEdge(pCFGNode& from, pCFGNode& to, string color, string style)
{
    // edgeouts << to.getCurAncsNode(curPSet).id() << "->" << to.id() << endl;
    edgeouts << from.id() << "->" << to.id() << "[style=" << style << "," << "color=" << color << "];" << endl;
}

void pCFGIterator::printNode(pCFGNode& pcfgn)
{
    //nodeouts << pcfgn.toDot() << endl;
    nodeouts << pcfgn.toDot();
}

void pCFGIterator::printSubGraphInit(unsigned int curPSet)
{
    stringstream subgraph_name;
    subgraph_name << "pset_" << curPSet;
    string subg_name = subgraph_name.str();
    nodeouts << "subgraph " << subg_name << "{\n";
    if(curPSet >= pSetColors.size()) {
        pSetColors.push_back(genRandColor(curPSet));
    }
    ROSE_ASSERT(curPSet >= 0 && curPSet < pSetColors.size());
    nodeouts << "node [color=" << pSetColors[curPSet] << "];\n";
}

void pCFGIterator::printSubGraphEnd()
{
    nodeouts << "}\n";
}

void pCFGIterator::writeToDot(string filename)
{
    ofstream file;
    file.open(filename.c_str());
    file << "digraph pcfg {\n";
    file << "node [shape = box]; " << endl;
    file << nodeouts.str() << edgeouts.str();
    file << "}\n";
    file.close();
}

string pCFGIterator::genRandColor(unsigned int curPSet)
{
    stringstream color;
    boost::mt19937 seed(static_cast<unsigned int> (time(0) + curPSet));
    boost::uniform_smallint<> color_component(0, 255);
    boost::variate_generator< boost::mt19937, boost::uniform_smallint<> > 
        hex(seed, color_component);
    int r = hex();
    int g = hex();
    int b = hex();
    color << "\"# " << HEX(r) << HEX(g) << HEX(b) << "\"";
    //cout << std::hex << r << endl;
    return color.str();
}
