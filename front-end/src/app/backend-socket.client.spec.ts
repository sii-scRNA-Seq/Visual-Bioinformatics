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
});
