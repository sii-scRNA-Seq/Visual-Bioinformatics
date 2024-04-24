import * as socketio from 'socket.io-client';
import { TestBed } from '@angular/core/testing';

import { BackendSocketClient } from './backend-socket.client';

describe('BackendSocketService', () => {
  let service: BackendSocketClient;

  beforeEach(() => {
    TestBed.configureTestingModule({});
    service = TestBed.inject(BackendSocketClient);
  });

  it('should be created', () => {
    expect(service).toBeTruthy();
  });

  describe('sendRequest', () => {
    it('should appropriately emit the given message', () => {
      const spy = jasmine.createSpyObj('Socket', ['connect', 'emit', 'on']);
      const ioSpy = jasmine.createSpy('io').and.returnValue(spy);
      spyOnProperty(socketio, 'io', 'get').and.returnValue(ioSpy);

      const message = {
        user_id: 'mock_user_id',
        blocks: [{ block_id: 'loaddata', foo: 42}],
      };
      service.sendRequest(message);
      expect(spy.emit).toHaveBeenCalled();
      //expect(spy).toHaveBeenCalledOnceWith("json", message);
    });
  });
});
