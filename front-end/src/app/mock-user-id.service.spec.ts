import { TestBed } from '@angular/core/testing';

import { MockUserIdService } from './mock-user-id.service';

describe('MockUserIdService', () => {
  let service: MockUserIdService;

  beforeEach(() => {
    TestBed.configureTestingModule({});
    service = TestBed.inject(MockUserIdService);
  });

  it('should be created', () => {
    expect(service).toBeTruthy();
  });
});
