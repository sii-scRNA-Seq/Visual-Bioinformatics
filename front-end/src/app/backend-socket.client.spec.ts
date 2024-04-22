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
      const message = {
        user_id: 'mock_user_id',
        blocks: [{ block_id: 'loaddata', foo: 42}],
      };
      const spy = spyOn(service, 'sendRequest');
      service.sendRequest(message);

      expect(true).toBe(false);
    });

    it('should break', () => {
      expect(true).toBe(false);
    });
  });
});
