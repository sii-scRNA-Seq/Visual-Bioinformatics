import { EventsMap, ReservedOrUserListener } from '@socket.io/component-emitter';
import { Socket } from 'socket.io-client';
import { TestBed } from '@angular/core/testing';

import { BackendSocketClient } from './backend-socket.client';
import { MockSocket } from './mock-socket';
import { SOCKET } from './socket';

describe('BackendSocketService', () => {
  let service: BackendSocketClient;

  beforeEach(() => {
    TestBed.configureTestingModule({
      providers: [
        { provide: SOCKET, useClass: MockSocket },
      ]
    });
  });

  it('should be created', () => {
    service = TestBed.inject(BackendSocketClient);
    expect(service).toBeTruthy();
  });

  describe('constructor', () => {
    it('should call socket.connect()', () => {
      const socket: Socket = TestBed.inject(SOCKET);
      const spy = spyOn(socket, 'connect');
      service = TestBed.inject(BackendSocketClient);
      expect(spy).toHaveBeenCalledOnceWith();
    });
  });

  describe('sendRequest', () => {
    beforeEach(()=>{
      service = TestBed.inject(BackendSocketClient);
    });

    it('should appropriately emit the given message', () => {
      const message = {
        user_id: 'mock_user_id',
        blocks: [{ block_id: 'loaddata', foo: 42}],
      };
      const socket: Socket = TestBed.inject(SOCKET);
      const spy = spyOn(socket, 'emit');
      service.sendRequest(message);
      expect(spy).toHaveBeenCalledOnceWith('json', message);
    });
  });

  describe('listen', () => {
    beforeEach(()=>{
      service = TestBed.inject(BackendSocketClient);
    });

    it('should make an appropriate call to socket.on()', () => {
      const foo = jasmine.createSpy('foo', (res: string) => { }); // eslint-disable-line @typescript-eslint/no-unused-vars
      const socket: Socket = TestBed.inject(SOCKET);
      const spy = spyOn(socket, 'on').and.callFake((str: string, callback: ReservedOrUserListener<EventsMap, EventsMap, string>) => {
        callback('test string!');
        return socket;
      });
      service.listen(foo);
      expect(spy).toHaveBeenCalledOnceWith('json', jasmine.any(Function));
      expect(foo).toHaveBeenCalledOnceWith('test string!');
    });
  });
});
