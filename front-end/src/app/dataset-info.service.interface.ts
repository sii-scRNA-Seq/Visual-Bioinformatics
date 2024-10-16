import { Observable } from 'rxjs';

import { DatasetInfo } from './dataset-info';

export interface DatasetInfoServiceInterface {
  datasetInfo: Observable<DatasetInfo[]>;

  setDatasetInfo(): Promise<void>;
}
