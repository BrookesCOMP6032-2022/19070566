import math
from dataclasses import dataclass

import networld


# a data container object for the taxi's internal list of fares. This
# tells the taxi what fares are available to what destinations at
# what price, and whether they have been bid upon or allocated. The
# origin is notably missing: that's because the Taxi will keep this
# in a dictionary indexed by fare origin, so we don't need to duplicate that
# here.

# Taxi fare and total taxi on service counter.
@dataclass
class TaxiGroupInformation:
    taxiCounter: int
    total: int
    taxisJobs: []
    avg: int
    daysReturns: []
    openFares: []
    time: int
    fareCount: int
    trafficProbabilityMap: []
    timeCost: []


group = TaxiGroupInformation(0, 0, [], 0, [], [], 0, 0, [], [])


class FareInfo:
    def __init__(self, destination, price):
        self.destination = destination
        self.price = price
        # bid is a ternary value: -1 = no, 0 = undecided, 1 = yes indicating whether this
        # taxi has bid for this fare.
        self.bid = 0
        self.allocated = False


''' A Taxi is an agent that can move about the world it's in and collect fares. All taxis have a
    number that identifies them uniquely to the dispatcher. Taxis have a certain amount of 'idle'
    time they're willing to absorb, but beyond that, they go off duty since it seems to be a waste
    of time to stick around. Eventually, they might come back on duty, but it usually won't be for
    several hours. A Taxi also expects a map of the service area which forms part of its knowledge
    base. Taxis start from some fixed location in the world. Note that they can't just 'appear' there:
    any real Node in the world may have traffic (or other taxis!) there, and if its start node is
    unavailable, the taxi won't enter the world until it is. Taxis collect revenue for fares, and 
    each minute of active time, whether driving, idle, or conducting a fare, likewise costs them £1.
'''


class Taxi:
    # message type constants
    FARE_ADVICE = 1
    FARE_ALLOC = 2
    FARE_PAY = 3
    FARE_CANCEL = 4

    '''constructor. The only required arguments are the world the taxi operates in and the taxi's number.
       optional arguments are:
       idle_loss - how much cost the taxi is prepared to absorb before going off duty. 256 gives about 4
       hours of life given nothing happening. Loss is cumulative, so if a taxi was idle for 120 minutes,
       conducted a fare over 20 minutes for a net gain to it of 40, then went idle for another 120 minutes,
       it would have lost 200, leaving it with only £56 to be able to absorb before going off-duty.
       max_wait - this is a heuristic the taxi can use to decide whether a fare is even worth bidding on.
       It is an estimate of how many minutes, on average, a fare is likely to wait to be collected.
       on_duty_time - this says at what time the taxi comes on duty (default is 0, at the start of the 
       simulation)
       off_duty_time - this gives the number of minutes the taxi will wait before returning to duty if
       it goes off (default is 0, never return)
       service_area - the world can optionally populate the taxi's map at creation time.
       start_point - this gives the location in the network where the taxi will start. It should be an (x,y) tuple.
       default is None which means the world will randomly place the Taxi somewhere on the edge of the service area.
    '''

    def __init__(self, world, taxi_num, idle_loss=256, max_wait=50, on_duty_time=0, off_duty_time=0, service_area=None,
                 start_point=None):
        self.total = 0
        self._world = world
        self.number = taxi_num
        self.onDuty = False
        self._onDutyTime = on_duty_time
        self._offDutyTime = off_duty_time
        self._onDutyPos = start_point
        self._dailyLoss = idle_loss
        self._maxFareWait = max_wait
        self._account = 0
        self._loc = None
        self._direction = -1
        self._nextLoc = None
        self._nextDirection = -1
        # this contains a Fare (object) that the taxi has picked up. You use the functions pickupFare()
        # and dropOffFare() in a given Node to collect and deliver a fare
        self._passenger = None
        # the map is a dictionary of nodes indexed by (x,y) pair. Each entry is a dictionary of (x,y) nodes that indexes a
        # direction and distance. such a structure allows rapid lookups of any node from any other.
        self._map = service_area
        if self._map is None:
            self._map = self._world.exportMap()
        # path is a list of nodes to be traversed on the way from one point to another. The list is
        # in order of traversal, and does NOT have to include every node passed through, if these
        # are incidental (i.e. involve no turns or stops or any other remarkable feature)
        self._path = []
        # pick the first available entry point starting from the top left corner if we don't have a
        # preferred choice when coming on duty
        if self._onDutyPos is None:
            x = 0
            y = 0
            while (x, y) not in self._map and x < self._world.xSize:
                y += 1
                if y >= self._world.ySize:
                    y = 0
                    x += self._world.xSize - 1
            if x >= self._world.xSize:
                raise ValueError("This taxi's world has a map which is a closed loop: no way in!")
            self._onDutyPos = (x, y)
        # this dict maintains which fares the Dispatcher has broadcast as available. After a certain
        # period of time, fares should be removed  given that the dispatcher doesn't inform the taxis
        # explicitly that their bid has not been successful. The dictionary is indexed by
        # a tuple of (time, originx, originy) to be unique, and the expiry can be implemented using a heap queue
        # for priority management. You would do this by initialising a self._fareQ object as:
        # self._fareQ = heapq.heapify(self._fares.keys()) (once you have some fares to consider)

        # the dictionary items, meanwhile, contain a FareInfo object with the price, the destination, and whether
        # or not this taxi has been allocated the fare (and thus should proceed to collect them ASAP from the origin)
        self._availableFares = {}

    # This property allows the dispatcher to query the taxi's location directly. It's like having a GPS transponder
    # in each taxi.
    @property
    def currentLocation(self):
        if self._loc is None:
            return -1, -1
        return self._loc.index

    # ___________________________________________________________________________________________________________________________
    # methods to populate the taxi's knowledge base

    # get a map if none was provided at the outset
    def importMap(self, newMap):
        # a fresh map can just be inserted
        if self._map is None:
            self._map = newMap
        # but importing a new map where one exists implies adding to the
        # existing one. (Check that this puts in the right values!)
        else:
            for node in newMap.items():
                neighbours = [(neighbour[1][0], neighbour[0][0], neighbour[0][1]) for neighbour in node[1].items()]
                self.addMapNode(node[0], neighbours)

                # incrementally add to the map. This can be useful if, e.g. the world itself has a set of

    # nodes incrementally added. It can then call this function on the existing taxis to update
    # their maps.
    def addMapNode(self, coords, neighbours):
        if self._world is None:
            return AttributeError("This Taxi does not exist in any world")
        node = self._world.getNode(coords[0], coords[1])
        if node is None:
            return KeyError("No such node: {0} in this Taxi's service area".format(coords))
        # build up the neighbour dictionary incrementally so we can check for invalid nodes.
        neighbourDict = {}
        for neighbour in neighbours:
            neighbourCoords = (neighbour[1], neighbour[2])
            neighbourNode = self._world.getNode(neighbour[1], neighbour[2])
            if neighbourNode is None:
                return KeyError(
                    "Node {0} expects neighbour {1} which is not in this Taxi's service area".format(coords, neighbour))
            neighbourDict[neighbourCoords] = (neighbour[0], self._world.distance2Node(node, neighbourNode))
        self._map[coords] = neighbourDict

    # ---------------------------------------------------------------------------------------------------------------------------
    # automated methods to handle the taxi's interaction with the world. You should not need to change these.

    # comeOnDuty is called whenever the taxi is off duty to bring it into service if desired. Since the relevant
    # behaviour is completely controlled by the _account, _onDutyTime, _offDutyTime and _onDutyPos properties of
    # the Taxi, you should not need to modify this: all functionality can be achieved in clockTick by changing
    # the desired properties.
    def comeOnDuty(self, time=0):
        if self._world is None:
            return AttributeError("This Taxi does not exist in any world")
        if self._offDutyTime == 0 or (self._onDutyTime <= time < self._offDutyTime):
            if self._account <= 0:
                self._account = self._dailyLoss
            self.onDuty = True
            onDutyPose = self._world.addTaxi(self, self._onDutyPos)
            group.taxiCounter += 1

            print("New taxi has just come on duty.\nThe total amount of taxis on duty now is: ", group.taxiCounter)
            self._nextLoc = onDutyPose[0]
            self._nextDirection = onDutyPose[1]

    # clockTick should handle all the non-driving behaviour, turn selection, stopping, etc. Drive automatically
    # stops once it reaches its next location so that if continuing on is desired, clockTick has to select
    # that action explicitly. This can be done using the turn and continueThrough methods of the node. Taxis
    # can collect fares using pickupFare, drop them off using dropoffFare, bid for fares issued by the Dispatcher
    # using transmitFareBid, and any other internal activity seen as potentially useful.
    def clockTick(self, world):
        # automatically go off duty if we have absorbed as much loss as we can in a day
        if self._account <= 0 and self._passenger is None:
            print("Taxi {0} is going off-duty".format(self.number))
            self.onDuty = False
            group.taxiCounter -= 1
            print("Taxi has gone off duty.\nThe current amount of taxi's on duty is:", group.taxiCounter)
            self._offDutyTime = self._world.simTime
        # have we reached our last known destination? Decide what to do now.
        if len(self._path) == 0:
            # obviously, if we have a fare aboard, we expect to have reached their destination,
            # so drop them off.
            if self._passenger is not None:
                if self._loc.dropoffFare(self._passenger, self._direction):
                    self._passenger = None
                # failure to drop off means probably we're not at the destination. But check
                # anyway, and replan if this is the case.
                elif self._passenger.destination != self._loc.index:
                    self._path = self._planPath(self._loc.index, self._passenger.destination)

        # decide what to do about available fares. This can be done whenever, but should be done
        # after we have dropped off fares so that they don't complicate decisions.
        faresToRemove = []
        for fare in self._availableFares.items():
            # remember that availableFares is a dict indexed by (time, originx, originy). A location,
            # meanwhile, is an (x, y) tuple. So fare[0][0] is the time the fare called, fare[0][1]
            # is the fare's originx, and fare[0][2] is the fare's originy, which we can use to
            # build the location tuple.
            origin = (fare[0][1], fare[0][2])
            # much more intelligent things could be done here. This simply naively takes the first
            # allocated fare we have and plans a basic path to get us from where we are to where
            # they are waiting.
            if len(self._path) == 0 and fare[1].allocated and self._passenger is None:
                # at the collection point for our next passenger?
                if self._loc.index[0] == origin[0] and self._loc.index[1] == origin[1]:
                    self._passenger = self._loc.pickupFare(self._direction)
                    # if a fare was collected, we can start to drive to their destination. If they
                    # were not collected, that probably means the fare abandoned.
                    if self._passenger is not None:
                        self._path = self._planPath(self._loc.index, self._passenger.destination)
                    faresToRemove.append(fare[0])
                # not at collection point, so determine how to get there
                else:
                    self._path = self._planPath(self._loc.index, origin)
            # get rid of any unallocated fares that are too stale to be likely customers
            elif self._world.simTime - fare[0][0] > self._maxFareWait:
                faresToRemove.append(fare[0])
            # may want to bid on available fares. This could be done at any point here, it
            # doesn't need to be a particularly early or late decision amongst the things to do.
            elif fare[1].bid == 0:
                if self._bidOnFare(fare[0][0], origin, fare[1].destination, fare[1].price):
                    self._world.transmitFareBid(origin, self)
                    fare[1].bid = 1
                else:
                    fare[1].bid = -1
        for expired in faresToRemove:
            del self._availableFares[expired]
        # may want to do something active whilst enroute - this simple default version does
        # nothing, but that is probably not particularly 'intelligent' behaviour.
        else:
            pass

        # the last thing to do is decrement the account - the fixed 'time penalty'. This is always done at
        # the end so that the last possible time tick isn't wasted e.g. if that was just enough time to
        # drop off a fare.
        self._account -= 1

    # called automatically by the taxi's world to update its position. If the taxi has indicated a
    # turn or that it is going straight through (i.e., it's not stopping here), drive will
    # move the taxi on to the next Node once it gets the green light.
    def drive(self, newPose):
        # if we have no direction yet, we need to get one
        if newPose[0] is None and newPose[1] == -1:
            self._nextLoc = None
            self._nextDirection = -1
        # if we do have a direction, as long as we are not stopping here,
        if self._nextLoc is not None:
            # and we have the green light to proceed,
            if newPose[0] == self._nextLoc and newPose[1] == self._nextDirection:
                nextPose = (None, -1)
                # vacate our old position and occupy our new node.
                if self._loc is None:
                    nextPose = newPose[0].occupy(newPose[1], self)
                else:
                    nextPose = self._loc.vacate(self._direction, newPose[1])
                # managed to occupy the destination location?
                if nextPose[0] == newPose[0] and nextPose[1] == newPose[1]:
                    self._loc = self._nextLoc
                    self._direction = self._nextDirection
                    # self._nextLoc = None
                    # self._nextDirection = -1
                    # not yet at the destination?
                    if len(self._path) > 0:
                        #  if we have reached the next path waypoint, pop it.
                        if self._path[0][0] == self._loc.index[0] and self._path[0][1] == self._loc.index[1]:
                            self._path.pop(0)
                        # otherwise continue straight ahead (as long as this is feasible)
                        else:
                            nextPose = self._loc.continueThrough(self._direction)
                            self._nextLoc = nextPose[0]
                            self._nextDirection = nextPose[1]
                            return
                # clear out the destination pose in preparation for computing a new one.
                self._nextLoc = None
                self._nextDirection = -1
        # Either get underway or move from an intermediate waypoint. Both of these could be
        # a change of direction. The logic above falls through to this test unless the 'straight ahead'
        # condition was met.
        if self._nextLoc is None and len(self._path) > 0:
            # we had better be in a known position!
            if self._loc.index not in self._map:
                raise IndexError("Fell of the edge of the world! Index ({0},{1}) is not in this taxi's map".format(
                    self._loc.index[0], self._loc.index[1]))
            #  if we are resuming from a previous path point, just pop the path
            if self._path[0][0] == self._loc.index[0] and self._path[0][1] == self._loc.index[1]:
                self._path.pop(0)
                # if after that, no more path is left, we're done. No need to drive further.
                if len(self._path) == 0:
                    return
            # we need to be going to a reachable location
            if self._path[0] not in self._map[self._loc.index]:
                raise IndexError(
                    "Can't get there from here! Map doesn't have a path to ({0},{1}) from ({2},{3})".format(
                        self._path[0][0], self._path[0][1], self._loc.index[0], self._loc.index[1]))
            # otherwise look up the next place to go from the map
            nextPose = self._loc.turn(self._direction, self._map[self._loc.index][self._path[0]][0])
            # update our next locations appropriately. If we're at the destination, or
            # can't move as expected, the next location will be None, meaning we will stop
            # here and will have to plan our path again.
            self._nextLoc = nextPose[0]
            self._nextDirection = nextPose[1]
            self._appendToTrafficProbability(self._path[0][0], self._path[0][1])

    # recvMsg handles various dispatcher messages.
    def recvMsg(self, msg, **args):
        timeOfMsg = self._world.simTime
        # A new fare has requested service: add it to the list of availables
        if msg == self.FARE_ADVICE:
            callTime = self._world.simTime
            self._availableFares[callTime, args['origin'][0], args['origin'][1]] = FareInfo(args['destination'],
                                                                                            args['price'])
            return
        # the dispatcher has approved our bid: mark the fare as ours
        elif msg == self.FARE_ALLOC:
            for fare in self._availableFares.items():
                if fare[0][1] == args['origin'][0] and fare[0][2] == args['origin'][1]:
                    if fare[1].destination[0] == args['destination'][0] and fare[1].destination[1] == \
                            args['destination'][1]:
                        fare[1].allocated = True
                        return

        # we just dropped off a fare and received payment, add it to the account
        elif msg == self.FARE_PAY:
            # Works out the actual time cost.
            # for route in group.timeCost:
            #     if route["destination"] == args["destination"]:
            #         route["withTraffic"] = abs(route["startTime"] - self._world.simTime) + route["noTraffic"]
            #         print(route["withTraffic"])

            self._account += args['amount']
            print("The current revenue: £", round(self._account, 2))
            return
        # a fare cancelled before being collected, remove it from the list
        elif msg == self.FARE_CANCEL:
            for fare in self._availableFares.items():
                if fare[0][1] == args['origin'][0] and fare[0][2] == args['origin'][1]:  # and fare[1].allocated:
                    del self._availableFares[fare[0]]
                    return

    # _____________________________________________________________________________________________________________________

    ''' HERE IS THE PART THAT YOU NEED TO MODIFY
    '''

    # this function should build your route and fill the _path list for each new
    # journey. Below is a naive depth-first search implementation. You should be able
    # to do much better than this!

    # This function uses the A* star algorithm to plan the path.
    def _planPath(self, origin, destination, **args):
        openNodes = {origin: [0, 0, 0, None]}  # Creates new open dic and adds origin value into dic
        closedNodes = {}
        #  While there is no undiscovered nodes.
        while len(openNodes) > 0:
            currentNodeF = None  # Sets the current nodes F value to a default of None.
            for node in openNodes:
                # If value has a lower f value set this value as current node.
                if currentNodeF is None or openNodes[node][0] < currentNodeF:
                    currentNode = node
                    # Update current F.
                    currentNodeF = openNodes[node][2]
                    # Stores the coordinates of the current node.
            currentCoordinates = openNodes[currentNode]
            # Removes from open end node list.
            del openNodes[currentNode]
            # Updates the list with the current lowest for ths node.
            closedNodes.update({currentNode: currentCoordinates[3]})
            if currentNode == destination:
                path = []
                parent = currentCoordinates[3]

                while parent is not None:
                    path.append(currentNode)
                    currentNode = parent
                    parent = closedNodes[currentNode]
                path.append(origin)
                path.reverse()
                return path
            else:
                # Updates all children nodes in the path.
                for child in (node for node in self._map[currentNode].keys() if node not in closedNodes):
                    # Calculates the cheapest path score.
                    childG = currentCoordinates[2] + self._world.distance2Node(
                        self._world.getNode(currentNode[0], currentNode[1]),
                        self._world.getNode(child[0], child[1]))

                    # Calculates the estimated heuristic score using Euclidean distance
                    # Problems with this method, nodes can get stuck
                    childH = (abs(destination[0] - child[0]) ** 2 + abs(destination[1] - child[1]) ** 2) * (
                            1 + self._getRoadProb(child))

                    # Current cheapest path.
                    childF = childG + childH
                    if child in openNodes:
                        if childG > openNodes[child][0]:
                            continue
                            # Updates  child node.
                    openNodes.update({child: [childG, childH, childF, currentNode]})
        return []

    #  This function appends values to a probabilistic map.
    def _appendToTrafficProbability(self, x, y):
        totalRoadUsers = sum([locI["trafficCount"] for locI in group.trafficProbabilityMap]) + 1
        flag = False
        for loc in group.trafficProbabilityMap:
            if loc["x"] == x and loc["y"] == y:
                # If in list increment traffic Counter
                loc["trafficCount"] += 1
                flag = True
                break

        if not flag:
            # If value is not found in probabilisticMap append new value.
            group.trafficProbabilityMap.append(
                {"x": x, "y": y, "probability": 1 / totalRoadUsers, "trafficCount": 1, "H": 0})

        for loc in group.trafficProbabilityMap:
            loc["probability"] = loc['trafficCount'] / totalRoadUsers

    # def _probabilisticPathPlanner(self, currentNode, previousNode):
    #     # Work out some costing value for a percentage
    #     for cross_section in group.trafficProbabilityMap:
    #         pass
    #     pass
    #
    #     # This method uses a Bayesian network of map crossroad probabilities, at every instance of a node travelling in a crossroad the map probabilities update.
    #
    #     # This function checks to see if the road already exists in the dictionary.
    #     def _check_if_road_exists_in_dictionary(self, keys, values):
    #         return all(
    #             [all([d[key] == value for key, value in zip(keys, values)]) for d in group.trafficProbabilityMap])
    #
    #     # This function gets the map probabilities of a specific road.
    #     # This will eventually be used in future uses where it can use it as a hash lookup table to find the shortest route using some sort of probability Heuristic

    # Gets the probability of a node using a dictionary search.
    def _getRoadProb(self, location):
        if len(group.trafficProbabilityMap) == 0:
            return 0
        for loc in group.trafficProbabilityMap:
            if loc["x"] == location[0] and loc["y"] == location:
                return loc["probability"]
        return 0

    def _bidOnFare(self, time, origin, destination, price):
        NoCurrentPassengers = self._passenger is None
        NoAllocatedFares = len([fare for fare in self._availableFares.values() if fare.allocated]) == 0
        TimeToOrigin = self._world.travelTime(self._loc, self._world.getNode(origin[0], origin[1]))
        TimeToDestination = self._world.travelTime(self._world.getNode(origin[0], origin[1]),
                                                   self._world.getNode(destination[1], destination[1]))

        group.timeCost.append({"noTraffic": TimeToDestination + TimeToOrigin,
                               "origin": origin,
                               "destination": destination,
                               "withTraffic": 0,
                               "startTime": self._world.simTime,
                               "probability": 0
                               })
        for fareTask in group.timeCost:
            if fareTask["withTraffic"] > 0:
                print("Time Cost:", (1 / fareTask["withTraffic"] - fareTask["noTraffic"]) / fareTask["probability"])
        # Check for other heuristic costs. and compare them again cost of probability * traffic probability.
        # probability = self._getCrossRoadProbabilty(origin)

        # heuristic + 1 * mean of all probability that occur in cross sections

        # timeTakenWithTraffic - timeWithoutTraffic * probability for all cross

        # withTraffic = noTraffic * 1 + probability

        FiniteTimeToOrigin = TimeToOrigin > 0
        FiniteTimeToDestination = TimeToDestination > 0
        CanAffordToDrive = self._account > TimeToOrigin
        FairPriceToDestination = price > TimeToDestination
        PriceBetterThanCost = FairPriceToDestination and FiniteTimeToDestination
        FareExpiryInFuture = self._maxFareWait > self._world.simTime - time
        EnoughTimeToReachFare = self._maxFareWait - self._world.simTime + time > TimeToOrigin

        SufficientDrivingTime = FiniteTimeToOrigin and EnoughTimeToReachFare
        WillArriveOnTime = FareExpiryInFuture and SufficientDrivingTime
        NotCurrentlyBooked = NoCurrentPassengers and NoAllocatedFares
        CloseEnough = CanAffordToDrive and WillArriveOnTime
        Worthwhile = PriceBetterThanCost and NotCurrentlyBooked
        Bid = CloseEnough and Worthwhile
        return Bid

# TODO

# this function decides whether to offer a bid for a fare. In general you can consider your current position, time,
# financial state, the collection and dropoff points, the time the fare called - or indeed any other variable that
# may seem relevant to decide whether to bid. The (crude) constraint-satisfaction method below is only intended as
# a hint that maybe some form of CSP solver with automated reasoning might be a good way of implementing this. But
# other metho dologies could work well. For best results you will almost certainly need to use probabilistic reasoning.
