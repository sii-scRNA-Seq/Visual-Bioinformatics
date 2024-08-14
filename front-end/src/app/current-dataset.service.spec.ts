import { TestBed } from '@angular/core/testing';

import { CurrentDatasetService } from './current-dataset.service';

describe('CurrentDatasetService', () => {
  let service: CurrentDatasetService;

  beforeEach(() => {
    TestBed.configureTestingModule({});
    service = TestBed.inject(CurrentDatasetService);
  });

  it('should be created', () => {
    expect(service).toBeTruthy();
  });
});
