import threading ,time

class without_import(threading.Thread):
    def run(self):
        time.sleep(.5)
        print 'worker: without done'

class with_import(threading.Thread):
    def run(self):
        import time
        time.sleep(.5)
        print 'worker: with done'

def test():
    ############# this works fine
    wo = without_import()
    wo.start()
    wo.join()
    print 'main: without import done'

    ######## without join(), thread also terminates normally
    w = with_import()
    w.start()
    print 'main: with import done'

    ##### this blocks indefinitely
    w = with_import()
    w.start()
    w.join()
    print 'main: with import done'

    print 'ALL TESTS PASS'

test()

