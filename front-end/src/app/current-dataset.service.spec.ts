import { TestBed } from '@angular/core/testing';

import { CurrentDatasetService } from './current-dataset.service';
import { DatasetInfoService } from './dataset-info.service';
import { MockDatasetInfoService } from './mock-dataset-info.service';

describe('CurrentDatasetService', () => {
  let service: CurrentDatasetService;
  // let datasetInfoService: DatasetInfoService;

  beforeEach(() => {
    TestBed.configureTestingModule({
      providers: [
        { provide: DatasetInfoService, useClass: MockDatasetInfoService },
      ]
    });
    service = TestBed.inject(CurrentDatasetService);
    // datasetInfoService = TestBed.inject(DatasetInfoService);
  });

  it('should be created', () => {
    expect(service).toBeTruthy();
  });
});
